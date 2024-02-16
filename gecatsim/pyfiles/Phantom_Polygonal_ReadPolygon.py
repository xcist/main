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


def extract_lesions(file_path):
    lesions = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        index = 0
        while index < len(lines):
            if index + 4 >= len(lines):
                break
            name = lines[index + 1].strip()
            material_index_str = lines[index + 2].strip()
            # try:
            #     material_index_parts = list(map(int, material_index_str.split()))
            # except ValueError:
            #     print(f"Error: Material index format '{material_index_str}' is not supported. ")
            #     return []
            num_polygons = int(lines[index + 3])
            line_vertices = []
            for i in range(num_polygons):
                polygon_line = lines[index + 4 + i].strip()
                vertex = tuple(map(str, polygon_line.split()))
                line_vertices.append(vertex)
            try:
                if len(line_vertices) != num_polygons:
                    raise ValueError(f"Error: Number of polygons doesn't match the stated count in line {index +3}")
            except ValueError as e:
                print(f"Error: {e}")
                return None
            lesions.append({
                'name': name,
                'material_index': material_index_str,
                'all_polygons' : line_vertices
            })
            # move to next lesion
            index += 4 + num_polygons
    return lesions

def extract_polygons(lesions):
    total_polygon_count = 0
    all_polygons_combined = []
    for lesion in lesions:
        each_lesion_polygons = lesion['all_polygons']
        all_polygons_combined.extend(each_lesion_polygons)
        total_polygon_count += len(each_lesion_polygons)

        for i, polygon in enumerate(all_polygons_combined, start=1):
            formed_polygon = " ".join([f"({point})" for point in polygon])
            # print(f"Polygon {i}: ({formed_polygon})")
    print(f"\nTotal polygons : {total_polygon_count}")



file_path ="female_10yr_lung_lesions.nrb"
lesions_ = extract_lesions(file_path)
polygons = extract_polygons(lesions_)
