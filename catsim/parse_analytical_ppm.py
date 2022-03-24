def parse_analytical_ppm(ppmPhantomFilename):
    #TODO: consider change obj to a class
    obj = dict.fromkeys(['materialList', 'center', 'half_axes', 'euler_angs', 'density', 'type', 'material', 'axial_lims', 'shape', 'clip'])
    with open(ppmPhantomFilename, 'r') as f:
        all_lines = f.readlines()

    key_indice = []
    for line in all_lines:
        if "materialList" in line:
            material_list = line.strip().rstrip(";").split('=')[-1].lstrip('{').rstrip('}').split()
            material_list = [tmp.strip().lstrip('{').rstrip('}').strip("'") for tmp in material_list]
            obj['materialList'] = material_list
        elif len(line) < 4: continue
        else:
            key_str, value_str = line.rstrip(";").split('=')
            key_idx = int(''.join(filter(str.isdigit, key_str)))
            key_indice.append(key_idx)

    key_indice = set(key_indice)
    for key in obj.keys():
        #if key != 'materialList': obj[key] = [[]]*len(key_indice)
        if key != 'materialList': obj[key] = [[] for _ in range(len(key_indice))]

    for line in all_lines:
        if len(line) < 4: continue
        elif "materialList" not in line:
            key_str, value_str = line.rstrip(";").split('=')
            key_word = key_str.split('.')[-1].split('(')[0].split('{')[0]
            key_word = str(key_word)
            #key_idx = int(''.joint([x for x in key_str.split('.')[-1] if x.isdigit()]))
            key_idx = int(''.join(filter(str.isdigit, key_str)))
            value_str = value_str.strip().rstrip(';').rstrip(']').lstrip('[').split(' ')
            #if len(value_str)==1:
            #    try:
            #        value = eval(value_str)
            #    except: breakpoint()
            #else:
            # TODO: currently all values will be saved as list, even for scalars. try a better way
            try: value = [eval(x) for x in value_str]
            except: value = []

            obj[key_word][key_idx-1] = value

    return obj

if __name__=="__main__":
    filename = "../phantom/FB_head.ppm"
    phobj = parse_analytical_ppm(filename)
    breakpoint()
