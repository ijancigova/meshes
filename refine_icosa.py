import numpy as np


def distance(a, b):
    return np.linalg.norm(np.array(a) - np.array(b))


def project_points(points):
    new_points = []
    for point in points:
        d = distance(point, [0.0, 0.0, 0.0])
        if d > epsilon:
            new_points.append(point / d)
        else:
            new_points.append(point)
    return new_points


def check_duplicate(point, point_id):
    for i, node in enumerate(nodes):
        if distance(node, point) < epsilon:
            return i
    return point_id


def refine(k, triangle):
    a = np.array(nodes[triangle[0]])
    b = np.array(nodes[triangle[1]])
    c = np.array(nodes[triangle[2]])
    ab = b - a
    ac = c - a
    bc = c - b
    point_id_list = []
    a_id = check_duplicate(a, len(nodes))
    point_id_list.append(a_id)
    row = k
    direction = 1   # direction=1 means to the right (bc direction), direction=-1 means to the left (ac direction)
    start = a
    while row > 0:
        for i in range(row):
            new = start + direction * ab / (k*1.0)
            new_id = check_duplicate(new, len(nodes))
            if new_id == len(nodes):               # only save new points
                nodes.append(list(new))
            point_id_list.append(new_id)
            start = new
        row -= 1
        if direction == 1:
            new = start + bc / (k*1.0)
        else:
            new = start + ac / (k*1.0)

        new_id = check_duplicate(new, len(nodes))
        point_id_list.append(new_id)
        if new_id == len(nodes):               # only save new points
            nodes.append(list(new))
        direction *= -1
        start = new

    # give the triangles their correct IDs
    triangles_template = [[int(v) - 1 for v in line.split()] for line in open('TriOrig'+str(k*k)+'.dat')]
    triangles = []
    for i in triangles_template:
        triangles.append([point_id_list[i[0]], point_id_list[i[1]], point_id_list[i[2]]])

    return triangles


t = (1+np.sqrt(5))/2.0
s = np.sqrt(1+t**2)
epsilon = 0.000001

nodes = [[t/s, 1/s, 0.0],
         [-t/s, 1/s, 0.0],
         [t/s, -1/s, 0.0],
         [-t/s, -1/s, 0.0],
         [1/s, 0.0, t/s],
         [1/s, 0.0, -t/s],
         [-1/s, 0.0, t/s],
         [-1/s, 0.0, -t/s],
         [0.0, t/s, 1/s],
         [0.0, -t/s, 1/s],
         [0.0, t/s, -1/s],
         [0.0, -t/s, -1/s]]

triangles = [[0, 8, 4], [0, 5, 10], [2, 4, 9], [2, 11, 5],
             [1, 6, 8], [1, 10, 7], [3, 9, 6], [3, 7, 11],
             [0, 10, 8], [1, 8, 10], [2, 9, 11], [3, 11, 9],
             [4, 2, 0], [5, 0, 2], [6, 1, 3], [7, 3, 1],
             [8, 6, 4], [9, 4, 6], [10, 5, 7], [11, 7, 5]]

next_triangles = []
for triangle in triangles:
    new_triangles = refine(2, triangle)
    next_triangles += new_triangles
triangles = next_triangles
nodes = project_points(nodes)

next_triangles = []
for triangle in triangles:
    new_triangles = refine(5, triangle)
    next_triangles += new_triangles
triangles = next_triangles
nodes = project_points(nodes)

# output mesh
output_file = open("sphere" + str(len(nodes)) + str("nodes.dat"), "w")
for p in nodes:
    output_file.write(str(p[0]) + " " + str(p[1]) + " " + str(p[2]) + "\n")
output_file.close()
output_file = open("sphere" + str(len(nodes)) + str("triangles.dat"), "w")
for t in triangles:
    output_file.write(str(int(t[0])) + " " + str(int(t[1])) + " " + str(int(t[2])) + "\n")
output_file.close()

# output also vtk for visualisation
output_file = open("sphere" + str(len(nodes)) + ".vtk", "w")
output_file.write("# vtk DataFile Version 3.0\n")
output_file.write("Data\n")
output_file.write("ASCII\n")
output_file.write("DATASET POLYDATA\n")
output_file.write("POINTS " + str(len(nodes)) + " float\n")
for p in nodes:
    output_file.write(str(p[0]) + " " + str(p[1]) + " " + str(p[2]) + "\n")
output_file.write("TRIANGLE_STRIPS " + str(len(triangles)) + " " + str(4 * len(triangles)) + "\n")
for t in triangles:
    output_file.write("3 " + str(int(t[0])) + " " + str(int(t[1])) + " " + str(int(t[2])) + "\n")
output_file.close()

exit()
