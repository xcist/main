import re
import numpy as np

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

            if '[' in value_str:
                # extract substring in []
                value_str = re.findall(r'\[([^]]*)\]', value_str)[0]
                split_valstr = value_str.rstrip(';').split(';')
                if len(split_valstr) == 1:
                    this_str = split_valstr[0].strip().rstrip(";").split(' ')
                    if len(this_str)==1 and this_str[0] == '':
                        value = []
                    else:
                        value = [eval(x) for x in this_str]
                else:
                    value = []
                    for this_splvalstr in split_valstr:
                        _tmp = this_splvalstr.split(' ')
                        #breakpoint()
                        this_value = [eval(x) for x in _tmp]
                        value.append(this_value)
            else:
                this_str = value_str.strip().rstrip(";")
                value = eval(this_str)
            #value_str = value_str.strip().rstrip(';').rstrip(']').lstrip('[').split(' ')
            #value_str = value_str.strip().rstrip(';').rstrip(']').lstrip('[').split(' ')
            #breakpoint()
            #value_str = re.split(' |;', value_str.strip().rstrip(';').rstrip(']').lstrip('['))
            # if there is only 1d in clip, make it 2d
            if 'clip' in key_word and len(value)>0 and len(np.array(value).shape)==1: value=[value]
            #if len(value_str)==1:
            #    try:
            #        value = eval(value_str)
            #    except: breakpoint()
            #else:
            # TODO: currently all values will be saved as list, even for scalars. try a better way
            #try: value = [eval(x) for x in value_str]
            #except: value = []

            obj[key_word][key_idx-1] = value

    #breakpoint()
    return obj

if __name__=="__main__":
    filename = "/projects/catsim/XCIST/Analytical_projector/phantom/FB_hip.ppm"
    phobj = parse_analytical_ppm(filename)
    breakpoint()
