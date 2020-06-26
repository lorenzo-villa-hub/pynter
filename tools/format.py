



def format_composition(string,all_math=False):
    new_string = []
    if all_math:
        new_string.append('$')
    sp = list(string)
    for c in sp:
        if c.isdigit():
            if sp.index(c) != len(sp)-1:
                nextc = sp[sp.index(c)+1] 
                if nextc.isdigit():
                    newc = c + nextc
                    skip = True
                else:
                    newc = c
                    skip = False
            else:
                newc = c
                skip = False
            if not skip:
                if all_math:
                    new_string.append('_{%s}'%newc)
                else:
                    new_string.append('$_{%s}$'%newc)
        else:
            new_string.append(c)
    if all_math:
        new_string.append('$')       
    new_string = ''.join(new_string)
    
    return new_string