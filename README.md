# Interpolations

Problem: in the step to locate the point within the tessellation, it would be extremely slow when the point actually lies outside. In this case, when we have n triangles after Delaunay Triangulation and m points outside the tessellation, the complexity becomes O(mn)
