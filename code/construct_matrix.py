
#----------------------construct the matrix from wgcna_pre_vs_post_otuid.txt--------------#
def pre_vs_post_otuid(path = "wgcna_pre_vs_post_otuid.txt"):

    f = open(path, 'r')

    r = 9 #row number
    c = 8 #column number

    # create the r by n empty matrix M
    matrix = []
    for i in range(r):
        matrix += [[]*c]

    # fill the matrix, each entry M_i_j is a list of otus shared between cluster i
    # (from the post) and cluster j (from the pre)
    for i, row in enumerate(f):
        list_of_entries = row.split(';')
        for entry in list_of_entries[:-1]:
            list_of_otus = entry.split(',')
            if list_of_otus[0] == '':
                matrix[i].append([])
            else:
                matrix[i].append(list_of_otus)
    return matrix

#----------------------construct the matrix from wgcna_pre_vs_post_count.txt--------------#
def pre_vs_post_count(path = "wgcna_pre_vs_post_count.txt"):

    f = open(path, 'r')

    r = 9 #row number
    c = 8 #column number

    # create the r by n empty matrix M
    matrix = []
    for i in range(r):
        matrix += [[]*c]

    # fill the matrix, each entry M_i_j is the number of otus shared between cluster i
    #(from the post) and cluster j (from the pre)
    for i, row in enumerate(f):
        list_of_entries = row.split('\t')
        for entry in list_of_entries[:-1]:
                matrix[i].append(entry)
    return matrix

#-------------------------------main---------------------------------#

if __name__ == "__main__":

    otuid_matrix = pre_vs_post_otuid()
    count_matrix = pre_vs_post_count()

    # Example: print the number otus that are shared between post-3 and pre-0
    print(otuid_matrix[3][0])

    #Example: print the list of otus that are shared between post-3 and pre-0
    print(count_matrix[3][0])
