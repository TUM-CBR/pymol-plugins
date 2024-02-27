import numpy as np
from numpy import float64
from numpy.typing import NDArray
from pymol import cmd

from typing import Dict, Iterable, List, NamedTuple, Tuple

class SpaceOfCubes(NamedTuple):
    x_size: int
    y_size: int
    z_size: int
    delta: float64
    vertices: NDArray[float64]

    def reshape_3d(self) -> NDArray[float64]:
        return self.vertices.reshape(
            (self.x_size, self.y_size, self.z_size, 3)
        )
    
    def all_indexes(self) -> Iterable[Tuple[int, int, int]]:
        for i in range(self.x_size):
            for j in range(self.y_size):
                for k in range(self.z_size):
                    yield (i,j,k)

def get_cube_mesh(
    v1: NDArray[float64],
    v2: NDArray[float64],
    side: float
) -> SpaceOfCubes:

    if len(v1) != 3 or len(v2) != 3:
        raise ValueError("Arrays must have shape (3)")
    
    # Calculate dimensions of the larger cube
    dimensions = np.abs(v2 - v1)
    
    # Calculate the number of unit cubes along each dimension
    num_cubes = (np.floor(dimensions) / side).astype(int)
    
    # Total number of vertices for all unit cubes (8 vertices per cube)
    total_vertices = np.prod(num_cubes)
    
    # Pre-allocate array for storing vertices. Shape: (total_vertices, 3)
    vertices = np.zeros((total_vertices, 3))
    
    # Counter for filling the vertices array
    counter = 0
    for x in range(num_cubes[0]):
        for y in range(num_cubes[1]):
            for z in range(num_cubes[2]):
                vertices[counter] = v1 + side*np.array([x, y, z])
                counter += 1
    
    return SpaceOfCubes(
        x_size=num_cubes[0],
        y_size=num_cubes[1],
        z_size=num_cubes[2],
        delta=np.float64(side),
        vertices=vertices
    )

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
    max_distance: float,
    state: int = 1,
    mesh_spacing: float = 1
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

    xyzs = [vertices[:,0], vertices[:,1], vertices[:,2]]

    bottom = np.array([np.min(v) for v in xyzs])
    top = np.array([np.max(v) for v in xyzs])

    mesh = get_cube_mesh(bottom, top, mesh_spacing)

    groups = find_regions(mesh.vertices, vertices, max_distance, mesh_spacing)

    def is_close(atoms: List[NDArray[float64]], group: NDArray[float64]):

        for atom in atoms:
            min_norm = np.min(np.linalg.norm(group - atom, axis=1))

            if min_norm <= max_distance + 1:
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