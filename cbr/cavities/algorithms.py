import itertools
from igraph import Graph, VertexClustering
import numpy as np
from numpy.lib import recfunctions
from numpy import int64, float64
from numpy.typing import NDArray
from open3d.geometry import PointCloud, Octree, OctreeInternalPointNode, OctreeNode, OctreeNodeInfo, OctreeLeafNode, PointCloud, VoxelGrid
from open3d.utility import Vector3dVector
import pandas as pd
from pymol import cmd
from typing import cast

from typing import Dict, List, NamedTuple

K_BOX_ID = "box_id"
K_BOX_SIZE = "size"
K_BOX_DEPTH = "depth"
K_BOX_X = "x"
K_BOX_Y = "y"
K_BOX_Z = "z"
K_BOX_CX = "cx"
K_BOX_CY = "cy"
K_BOX_CZ = "cz"

class ProteinCavity(NamedTuple):
    boxes: pd.DataFrame

class FindCavitiesGraph(NamedTuple):
    graph: Graph
    depth: int

class FindCavitiesResult(NamedTuple):
    context: 'FindCavitiesContext'
    graphs: List[FindCavitiesGraph]
    nodes_df: pd.DataFrame
    edges_df: pd.DataFrame
    depths_to_volume: Dict[int, float]

    def __get_boxes(self, graph: Graph) -> ProteinCavity:
        node_ids = cast(List[int], [node['name'] for node in graph.vs])
        nodes_df = self.nodes_df.loc[node_ids]

        return ProteinCavity(nodes_df)

    def __get_cavities_for(self, cavities: FindCavitiesGraph, min_volume: int, max_volume: int):

        graph = cavities.graph
        unit_volume = self.depths_to_volume[cavities.depth]
        min_units = min_volume / unit_volume
        max_units = max_volume / unit_volume

        groups = cast(VertexClustering, graph.connected_components())
        accepted = [
            self.__get_boxes(graph.vs[group].subgraph())
            for group in groups
            if len(group) <= max_units and len(group) >= min_units
        ]

        return accepted

    def get_cavities(self, min_volume: int = 2, max_volume: int = 1000) -> List[ProteinCavity]:
        values = [
            cavity
            for graph in self.graphs
            for cavity in self.__get_cavities_for(graph, min_volume, max_volume)
        ]
        return values

class FindCavitiesContext(NamedTuple):
    points: NDArray[float64]
    point_cloud: PointCloud
    octree: Octree

    @staticmethod
    def construct(points: NDArray[float64]) -> 'FindCavitiesContext':

        pcd = PointCloud()
        pcd.points = Vector3dVector(points)

        # For protein structures, this leads to the smallest boxes
        # measuring almost 1A. This is sufficient for most cases.
        octree = Octree(max_depth=6)
        octree.convert_from_point_cloud(pcd)

        return FindCavitiesContext(
            points=points,
            point_cloud=pcd,
            octree=octree
        )

    def get_empty_corners(self) -> pd.DataFrame:
        """Constructs a DataFrame that contains the coordinates of the corners of all of the
        nodes of the Octree which are empty."""

        corners: List[NDArray[float64]] = []
        N_CORNERS = 8

        # Enumerate children coordinates just like
        # Open3d does the enumeration
        # https://github.com/isl-org/Open3D/blob/f5f672b4af1fc81e423c3c1b6215497f5a8816ea/cpp/open3d/geometry/Octree.cpp#L700
        children_coords_offset = np.array([
            [i % 2, int(i/2) % 2, int(i/4) % 2]
            for i in range(N_CORNERS)
        ])

        to_center = np.array([1,1,1])
        box_id = 0

        def apply(node: OctreeNode, node_info: OctreeNodeInfo):
            nonlocal corners, box_id

            child_size = node_info.size / 2

            if isinstance(node, OctreeInternalPointNode):
                for i,child in enumerate(node.children):
                    if child is not None:
                        continue

                    box_id += 1
                    uid = box_id
                    child_origin = node_info.origin + children_coords_offset[i]*child_size
                    center = child_origin + child_size/2 * to_center

                    # Convert corner's floating point coordinates as an int. We only need the first
                    # decimal position as that gives you a resolution of 0.1A, more than enough for our
                    # purposes. Representing the positions as integers allows one to quickly align corners
                    # using data manipulation libraries
                    box_vertices = ((np.tile(child_origin, (N_CORNERS, 1)) + children_coords_offset*child_size)*10).astype('int')

                    # We create an array with the fields common to every record, then we repeated
                    # once for every vertex of the box (8 times in 3d space)
                    common_values = np.repeat(
                        [
                            np.concatenate([
                                [uid, child_size, node_info.depth + 1],
                                center
                            ])
                        ],
                        len(box_vertices), axis=0
                    )

                    corners.append(np.hstack([common_values, box_vertices]))


        self.octree.traverse(apply)
        values = np.concatenate(corners)

        return pd.DataFrame({
            K_BOX_ID: values[:,0].astype('int'),
            K_BOX_SIZE: values[:,1],
            K_BOX_DEPTH: values[:,2].astype('int'),
            K_BOX_X: values[:,3],
            K_BOX_Y: values[:,4],
            K_BOX_Z: values[:,5],
            K_BOX_CX: values[:,6].astype('int'),
            K_BOX_CY: values[:,7].astype('int'),
            K_BOX_CZ: values[:,8].astype('int')
        })
    
    def construct_graph(self, corners: pd.DataFrame) -> FindCavitiesResult:

        CORNER_JOIN_KEYS = [K_BOX_CX, K_BOX_CY, K_BOX_CZ, K_BOX_DEPTH]
        K_BOX_ID_FROM = "box_id_from"
        K_BOX_ID_TO = "box_id_to"

        # Match the corners of each box with corners of other boxes
        # of the same depth. As all boxes come from an octree, adjecent
        # boxes will always share corners
        aligned = corners.rename({K_BOX_ID: K_BOX_ID_FROM}, axis=1) \
            .merge(
                corners.rename({K_BOX_ID: K_BOX_ID_TO}, axis=1),
                on=CORNER_JOIN_KEYS
            )
        
        # Filter out boxes matched to themselves
        K_COUNT_DUMMY = "count"
        aligned = aligned[aligned[K_BOX_ID_FROM] != aligned[K_BOX_ID_TO]]
        aligned[K_COUNT_DUMMY] = 1

        # Each edge will appear twice but in reverse order. We ensure that
        # the edge with the lower identifier always appears in the 'from'
        # and the edge with the higher identifier in the 'to'
        # This operation will create duplicated rows.
        min_of_edges = aligned[[K_BOX_ID_FROM, K_BOX_ID_TO]].min(axis=1)
        max_of_edges = aligned[[K_BOX_ID_FROM, K_BOX_ID_TO]].max(axis=1)
        aligned[K_BOX_ID_FROM] = min_of_edges
        aligned[K_BOX_ID_TO] = max_of_edges

        # Two boxes (of the same depth) are connected only if they have 4 corners in common.
        # However, we created duplicated rows in the previous step, so the actual magic value
        # is 8 = 2*4
        edges = aligned.groupby(by=[K_BOX_ID_FROM, K_BOX_ID_TO, K_BOX_DEPTH])[K_COUNT_DUMMY].count()
        edges = edges[edges == 8]
        edges_df = edges.index.to_frame()
        min_depth = cast(int, edges_df[K_BOX_DEPTH].min())
        max_depth = cast(int, edges_df[K_BOX_DEPTH].max())

        depths = (corners.groupby(K_BOX_DEPTH)[K_BOX_SIZE].max() ** 3).to_dict()

        nodes_df = corners.groupby(K_BOX_ID)[[K_BOX_X, K_BOX_Y, K_BOX_Z, K_BOX_SIZE]].max()

        graphs = [
            FindCavitiesGraph(
                Graph.DataFrame(
                    edges_df[edges_df[K_BOX_DEPTH] == d].drop(K_BOX_DEPTH, axis=1),
                    use_vids=False,
                    directed=False
                ),
                d
            )
            for d in range(min_depth, max_depth + 1)
        ]

        return FindCavitiesResult(
            self,
            graphs,
            nodes_df,
            edges_df,
            depths
        )

    @staticmethod
    def find_cavities(points: NDArray[float64]):
        ctx = FindCavitiesContext.construct(points)
        corners = ctx.get_empty_corners()
        return ctx.construct_graph(corners)

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