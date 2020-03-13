
###############################################################################


###############################################################################
def change_file (input_file , output_file=None, back_up_file = True,
                 str_to_modify = {}, lines_to_add=[], check_added_lines = True, lines_to_remove = []):
    
    ''' function to generate a copy of a file changing a string in specific lines
         In case input and output files are the same a copy of the old file is 
         created with name 'old_{input_file}'
    
     input_file: input file to change
     output_file: output file to generate
                  - Default is 'None' - in this case the same is set equal as the input file
                  - In case output_file and input_file have the same name a copy of 
                    the original input file is saved
     str_to_modify: dictionary of strings to be modified - format: {'old string':'new string'}
     lines_to_add: list of lines to be added at the end of the file: format: ['string 1','string 2']
     check_added_lines: - True - if line is already present is not added
                        - False - line is added indipendently of the others
    ''' 
    import re
    import os

    origin_file = open(input_file,'r')
    lines_list = origin_file.readlines()
    origin_file.close()
    
    if output_file == None or input_file == output_file:
        output_file = input_file
        if back_up_file == True:
            os.rename(input_file, 'old_' + input_file) 
        else:
            os.remove(input_file)
    

    new_file = open(output_file,'w')
         
      # searching lines in input file   
    for line in lines_list:
        
        rem_line = False
        mod_line = False 
    # checking lines to remove 
        for l in lines_to_remove:
            if l + '\n' == line:
                rem_line = True
        for string in str_to_modify: 
            # emulating 'grep' command
            target_line = re.findall(string, line)
            # if line is the target line replace string                      
            if target_line != [] and rem_line == False:
                # writing line with new string
                new_file.write(line.replace(string,str_to_modify[string]))
                mod_line = True                  
         # if no modified line has been written and line is not to be removed write original line
        if mod_line == False and rem_line == False:
             new_file.write(line)   
              # add lines 
         
    for l in lines_to_add:
          if check_added_lines:              
              if l + '\n' not in lines_list:            
                  new_file.write(l + '\n')
          else:
              new_file.write(l + '\n')
    
    new_file.close()
    
###############################################################################
           