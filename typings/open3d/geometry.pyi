from numpy import float64
from numpy.typing import NDArray
from typing import Any, Callable, Iterable, Optional
from open3d.utility import Vector3dVector

class PointCloud:

    def __init__(self) -> None: ...

    @property
    def points(self) -> Vector3dVector: ...

    @points.setter
    def points(self, points: Vector3dVector) -> None: ...

class OctreeNode:
    ...

class OctreeInternalPointNode(OctreeNode):

    @property
    def children(self) -> Iterable[Optional[OctreeNode]]: ...

class OctreeNodeInfo:

    @property
    def size(self) -> int: ...

    @property
    def origin(self) -> NDArray[float64]: ...

    @property
    def depth(self) -> int: ...

class Octree:

    def __init__(self, max_depth: int = ...) -> None: ...

    def convert_from_point_cloud(self, pcd: PointCloud) -> None: ...

    def traverse(self, fn: Callable[[OctreeNode, OctreeNodeInfo], Any]) -> None: ...