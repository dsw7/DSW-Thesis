# Written by dsw7@sfu.ca
# A function for generating custom homemade Avogadro compliant .xyz files


def get_XYZ_file(input_list, file_name):
    # I inherited this function from the now defunct PDBFileDriver library which
    # I wrote in early 2017 - still has many useful functions
    """
    PDBFileDriver.get_XYZ_file(input_list)
    
    A base script for generating .xyz Avogadro input files.
    
    Parameters : input_list, a list of elements and 3d coordinates formatted as follows:
                 input_list = [
                               ['N',1.00,2.00,1.00],
                               ['H',1.00,4.00,1.00], ...
                              ]
                 file_name, name of export file.
    
    The script generates:
        xyz_data.xyz
    """
    file = open(file_name, 'w')   
    file.write(str(len(input_list)) + '\n')
    file.write('PDBFileDriver generated XYZ file\n')    
    for i in input_list:
        string = ' {} {} {} {} \n'.format(*i)
        file.write(string)
    file.close()
    
    

