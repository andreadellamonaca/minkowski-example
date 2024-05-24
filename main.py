import minkowski

def read_points_from_file(file_path):
	points = []
	with open(file_path, 'r') as file:
		for line in file:
			x, y = map(int, line.strip().split(','))
			points.append((x, y))
	return points
	
a = read_points_from_file("./384_A.txt")
b = read_points_from_file("./384_B.txt")

q = minkowski.minkowski_sum(a,b)

print(f"N points: {len(q)}")
