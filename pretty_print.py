def print_list(l):
	for value in l:
		print("%0.2f" % value)

def print_dictionary(dic):
	for key, value in dic.items():
		if type(value) is list:
			l = [round(i,2) for i in value]
			print(key, l)
		elif type(value) is dict:
			l = [(k,round(v,2)) for k,v in value.items()]
			print(key, l)
		else:
			print("%s,%s") % (key,value)

def print_list_to_file(l, filename):
	fwriter = open(filename, 'w+')
	for value in l:
		fwriter.write("%s\n" % value)
	fwriter.close()

def print_dictionary_to_file(dic, filename, header=None):
	fwriter = open(filename, 'w')

	if header:
	    fwriter.write("%s\n" % header)

	for key, value in dic.items():
            if type(value) is list:
	        
                fwriter.write("%s" % key)
	        
                for item in value:
                    
	                fwriter.write(",%s" % item)
	                fwriter.write("\n")
                    
            else:
	            fwriter.write("%s,%s\n" % (key,value))

	fwriter.close()
