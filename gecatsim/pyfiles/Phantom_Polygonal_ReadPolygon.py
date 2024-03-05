from gecatsim.pyfiles.CommonTools import *

def Phantom_Polygonal_ReadPolygon(Verts):
    ddir = my_path.find_dir("phantom", "poly_bin")
    filename = 'poly{}'.format(Verts)
    with open(os.path.join(ddir, filename),'rb') as fid:
        data_array = np.fromfile(fid, dtype=np.int32, count=4)
        nv_sz1, nv_sz2, vx_sz1, vx_sz2 = data_array
        tmp = np.fromfile(fid, dtype=np.float64)
        nV = tmp[:nv_sz1*nv_sz2].reshape(nv_sz2,nv_sz1).T
        Vx = tmp[nv_sz1*nv_sz2:].reshape(vx_sz2,vx_sz1).T

    return Vx,nV

def extract_polygonal_objects(file_path):
    objects = []
    with open(file_path, 'r') as file:
        lines = [line.strip() for line in file.readlines()]
        index = 0
        combined_polygons = []
        while index < len(lines):
            if index + 4 >= len(lines):
                break
            object_name = lines[index + 1]
            material_index = int(lines[index + 2])
            try:
                num_polygons = int(lines[index + 3])
            except ValueError:
                print(f"Error: Number of polygons doesn't match the stated count in line {index + 3}")
                break

            all_polygon_vertices = []

            for i in range(num_polygons):
                read_line_of_coordinates = lines[index + 4 + i]
                each_polygon_coordinates = [tuple(map(float, point.split(','))) for point in read_line_of_coordinates.split()]
                all_polygon_vertices.append(each_polygon_coordinates)
            combined_polygons.extend(all_polygon_vertices)

            objects.append({
                'object_name': object_name,
                'material_index': material_index,
                'number_of_polygons': len(all_polygon_vertices),
                'polygons_per_object': all_polygon_vertices,
                # 'combined_polygons': combined_polygons    # If all the polygons shall be passed together
            })
            # move to next objects
            index += 4 + num_polygons

    return objects


file_path ="../phantom/reduced_female_10yr_lung_lesions.nrb"
_objects = extract_polygonal_objects(file_path)
print(f"\nObject names: ")
for obj in _objects:
    print(obj['object_name'])

