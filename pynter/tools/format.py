



def format_composition(string,all_math=False):
    
    def write(newc,new_string,skip,all_math):
        if not skip:
            if all_math:
                    new_string.append('_{%s}'%newc)
            else:
                new_string.append('$_{%s}$'%newc)
        return new_string
            
    skip = False
    new_string = []
    if all_math:
        new_string.append('$')
    sp = list(string)
    
    for c in sp:
        
        if c.isdigit():
            try:
                nextc = sp[sp.index(c)+1]
            except:
                nextc = ''
            if nextc.isdigit():
                newc = c + nextc
                new_string = write(newc,new_string,skip,all_math)
                skip = True
            else:
                newc = c
                new_string = write(newc,new_string,skip,all_math)
                skip = False

        else:
            new_string.append(c)
    if all_math:
        new_string.append('$')       
    new_string = ''.join(new_string)
    
    return new_string