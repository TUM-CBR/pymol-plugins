import itertools
import numpy as np
import pandas as pd
from numpy.lib import recfunctions
from numpy import int64, float64
from numpy.typing import NDArray
from open3d.geometry import Octree, OctreePointColorLeafNode, OctreeInternalNode, OctreeInternalPointNode, OctreeLeafNode, PointCloud, VoxelGrid
from open3d.utility import Vector3dVector
from pymol import cmd

from typing import Dict, List, NamedTuple

def get_corners(tree: Octree):

    xs = []
    ys = []
    zs = []
    nids = []
    depths = []
    sizes = []
    counts = []


    node_id = 0

    spacing = [0, 1]
    coord_items = [
        (x,y,z)
        for x in spacing
        for y in spacing
        for z in spacing
    ]

    coords = np.array(coord_items)

    def apply(node, node_info):
        nonlocal xs, ys, zs, nids, depths, node_id, sizes

        node_id += 1
        corners = coords * node_info.size + node_info.origin

        for corner in corners:
            origin = tuple(np.round(i,2) for i in corner)

            nids.append(node_id)
            depths.append(node_info.depth)
            sizes.append(node_info.size)
            xs.append(origin[0])
            ys.append(origin[1])
            zs.append(origin[2])

            if hasattr(node, 'indices'):
                counts.append(len(getattr(node, 'indices')))
            else:
                counts.append(0)

    tree.traverse(apply)

    return pd.DataFrame({
        "xs": xs,
        "ys": ys,
        "zs": zs,
        "nids": nids,
        "depths": depths,
        "sizes": sizes,
        "points": counts
    })


class IndexEntry(NamedTuple):
    value: float64
    i: int

class IndexedVertices(NamedTuple):
    sorted_axes: List[NDArray[float64]]
    sorted_to_vertex: List[NDArray[int64]]
    vertices: NDArray[float64]
    voxels: VoxelGrid

    @staticmethod
    def create(vertices: NDArray[float64]) -> 'IndexedVertices':
        (_, dims) = vertices.shape
        sorted_axes = [
            np.sort(vertices[:,d])
            for d in range(dims)
        ]
        sorted_to_vertex = [
            np.argsort(vertices[:,d])
            for d in range(dims)
        ]

        pcd = PointCloud()
        pcd.points = Vector3dVector(vertices)
        voxels = VoxelGrid.create_from_point_cloud(pcd, voxel_size=1)

        return IndexedVertices(
            sorted_axes=sorted_axes,
            sorted_to_vertex=sorted_to_vertex,
            vertices=vertices,
            voxels=voxels
        )
    
    @property
    def bottom(self) -> NDArray[float64]:
        return np.array([
            sorted_axis[0]
            for sorted_axis in self.sorted_axes
        ])
    
    @property
    def top(self) -> NDArray[float64]:
        return np.array([
            sorted_axis[-1]
            for sorted_axis in self.sorted_axes
        ])
    
    @property
    def dims(self) -> int:
        return len(self.sorted_axes)
    
    def empty_spaces_mesh(
        self,
        cutoff_distance: float,
        stride: float
    ) -> NDArray[float64]:
        
        sorted_axes = self.sorted_axes
        delta_vertex = np.array([cutoff_distance for _ in sorted_axes])

        bottom = self.bottom
        top = self.top

        # Calculate dimensions of the larger cube
        dimensions = np.abs(top - bottom)
        
        # Calculate the number of unit cubes along each dimension
        num_cubes = (np.floor(dimensions) / stride).astype(int)

        mesh = np.array([
            bottom + stride*np.array([x,y,z])
            for x in range(num_cubes[0])
            for y in range(num_cubes[1])
            for z in range(num_cubes[2])
        ])

        included_test = self.voxels.check_if_included(Vector3dVector(mesh))

        mesh_low = mesh - delta_vertex
        mesh_high = mesh + delta_vertex

        bottom_indices = [
            np.searchsorted(axis, mesh_low[:,d])
            for d,axis in enumerate(sorted_axes)
        ]

        top_indices = [
            np.searchsorted(axis, mesh_high[:,d], side='right')
            for d,axis in enumerate(sorted_axes)
        ]

        def mask_vertices(bottom: NDArray[int64], top: NDArray[int64]) -> bool:
            items = [
                self.sorted_to_vertex[d][bottom[d]:top[d]]
                for d in range(self.dims)
            ]

            base = items[0]

            for item in items[1:]:
                base = base[np.isin(base, item)]

            return len(base) == 0
        
        stride = 30000
        values = []
        for i in range(0, len(mesh), stride):
            top_i = i + stride
            indexes = [
                recfunctions.unstructured_to_structured(
                    np.concatenate([
                        np.stack(
                            [
                                row,
                                np.repeat(i, len(row))
                            ],
                            axis=1
                        )
                        for i,(bottom_ix, top_ix) in enumerate(zip(bottom_indices[d][i:top_i], top_indices[d][i:top_i]))
                        for row in [self.sorted_to_vertex[d][bottom_ix:top_ix]]
                    ]),
                    dtype=[('index', int64), ('mesh', int64)]
                )
                for d in range(self.dims)
            ]
            base = indexes[0]

            for item in indexes[1:]:
                base = base[np.isin(base, item)]

            values.append(base)
        
        mask_vectorized = np.vectorize(mask_vertices, signature=f"({self.dims}),({self.dims})->()")

        mask = mask_vectorized(
            np.stack(bottom_indices).transpose(),
            np.stack(top_indices).transpose()
        )
        return mesh[mask]

def find_nearest_distances(vertices: NDArray[float64], values: NDArray[float64]) -> NDArray[float64]:
    result = np.zeros(len(vertices)) + np.inf

    for value in values:
        candiates = np.stack([result, np.linalg.norm(vertices - value, axis=1)], axis=1)
        result = np.min(candiates, axis=1)

    return result

def make_groups(
    vertices: NDArray[float64],
    max_distance: float
) -> List[NDArray[float64]]:
    
    groups: List[NDArray[float64]] = []

    for i,vertex in enumerate(vertices):

        if len(groups) == 0:
            groups = [np.array([vertex])]
            continue

        try:
            candidates = [
                i
                for i,group in enumerate(groups)
                for closest_distance in [np.min(np.linalg.norm(group - vertex, axis=1))]
                    if closest_distance <= max_distance
            ]
        except Exception as e:
            raise e

        if len(candidates) > 1:
            # The current vertex joins two groups together
            # combine them into one

            # Candidates are indexes of groups in ascending order (due to the use of "enumerate")
            # We reverse them so we can pop them by index
            candidates.reverse()
            new_values = [groups.pop(i) for i in candidates] + [np.array([vertex])]
            new_group = np.concatenate(new_values, axis=0)
            groups.append(new_group)
        elif len(candidates) == 1:
            # vertex will be added to an existing group
            i = candidates[0]
            groups[i] = np.append(groups[i], [vertex], axis=0)
        else:
            # vertex belongs to a new group
            groups.append(np.array([vertex]))


    return groups

def find_regions(
    mesh: NDArray[float64],
    structure: NDArray[float64],
    max_distance: float,
    mesh_spacing: float
) -> List[NDArray[float64]]:
    distances = find_nearest_distances(mesh, structure)
    indexes = [i for i,d in enumerate(distances) if d > max_distance]
    return make_groups(
        mesh[indexes],
        mesh_spacing
    )

class Cavities(NamedTuple):
    selection: str
    groups: List[NDArray[float64]]
    groups_by_resv: List[List[int]]

def find_cavities(
    selection: str,
    cutoff_distance: float,
    state: int = 1,
    stride: float = 1
):
    vertices_dict: Dict[int, List[NDArray[float64]]] = {}

    def apply(resv:int, x: float, y: float, z: float):
        nonlocal vertices_dict

        vic = np.array([x,y,z])
        if resv not in vertices_dict:
            vertices_dict[resv] = [vic]
        else:
            vertices_dict[resv].append(vic)

    cmd.iterate_state(
        state,
        selection,
        'apply(resv, x, y, z)',
        space={'apply': apply}
    )

    vertices = np.array([
        vertex
        for group in vertices_dict.values()
        for vertex in group
    ])

    indexed_structure = IndexedVertices.create(vertices)

    mesh = indexed_structure.empty_spaces_mesh(cutoff_distance, stride)

    groups = find_regions(mesh, vertices, cutoff_distance, stride)

    def is_close(atoms: List[NDArray[float64]], group: NDArray[float64]):
        nonlocal cutoff_distance

        for atom in atoms:
            min_norm = np.min(np.linalg.norm(group - atom, axis=1))

            if min_norm <= cutoff_distance + 1:
                return True
        return False

    groups_by_resv = [
        [resv for resv, atoms in vertices_dict.items() if is_close(atoms, group)]
        for group in groups
    ]

    return Cavities(
        selection = selection,
        groups = groups,
        groups_by_resv = groups_by_resv
    )