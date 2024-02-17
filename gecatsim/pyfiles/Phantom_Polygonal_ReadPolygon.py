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

def extract_objects(file_path):
    objects = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        index = 0
        all_polygons = []
        while index < len(lines):
            if index + 4 >= len(lines):
                break
            object_name = lines[index + 1].strip()
            material_index_str = lines[index + 2].strip()
            material_index_parts = []
            try:
                material_index_parts = list(map(int, material_index_str.split()))
            except ValueError:
                print(f"Error: Material index format '{material_index_str}' is not supported. ")

            try:
                num_polygons = int(lines[index + 3])
            except ValueError:
                print(f"Error: Number of polygons doesn't match the stated count in line {index + 3}")
                break

            all_polygon_vertices = []

            for i in range(num_polygons):
                read_line_of_coordinates = lines[index + 4 + i].strip()
                each_polygon_coordinates = tuple(map(str, read_line_of_coordinates.split()))
                all_polygon_vertices.append(each_polygon_coordinates)
            all_polygons.extend(all_polygon_vertices)

            objects.append({
                'object_name': object_name,
                'material_index': material_index_parts,
                'polygons_per_object' : all_polygon_vertices
            })
            # move to next objects
            index += 4 + num_polygons

    return objects, all_polygons


file_path ="female_10yr_lung_lesions.nrb"
_objects, _all_polygons = extract_objects(file_path)
print(f"\n Total number of extracted polygons: {len(_all_polygons)}")

