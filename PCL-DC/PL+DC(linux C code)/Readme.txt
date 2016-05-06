IPHITS2.cpp

usage ./IPHITS2 cmd_line_file link_file content_file result_file

in cmd_line_file:
    Link link_file_name  %(link file name )
    Attribute cotent_file_name %( content file name)
    n number of nodes
    D number of attributes
    K number of commmunities
    ep precision of EM algorithm(0.000001)
    N number of iterations(1000)
    lamda value of lamda parameter(>=0)
    verbosity dispaly level(0,1,2,3)
    seed random seed(1)

link_file: the file stores the adjacent matrix
content_file: the file stores the attributes 
result_file: the fiel will store the returned result
