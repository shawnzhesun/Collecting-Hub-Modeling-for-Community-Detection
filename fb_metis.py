import networkx as nx
import sys

def convertGraph(graph, output):

	output.write(str(len(graph.node)) + " " + str(len(graph.edges())) + "\n")
	for i in range(1, len(graph.node) + 1):
		# print str(i) + "\t"
		if (str(i) in graph.nodes()):
			output.write(" ".join(graph.neighbors(str(i))) + "\n")
		else:
			output.write('\n')



def main():
	filename = sys.argv[1]+'.edges'
	print(filename)

	file = open(filename,'r')

	G = nx.Graph();

	for line in file:
		nodes = line.rstrip('\n').split(' ')
		# print nodes[0]+' '+nodes[1]
		G.add_edge(nodes[0],nodes[1])


	out_filename = sys.argv[1]+'.metis'
	out_file = open(out_filename,'w')

	convertGraph(G,out_file)

	feat_filename = sys.argv[1]+'.feat'
	feat_file = open(feat_filename)
	feat_out = open(sys.argv[1]+'.content','w')


	feat_count = len(feat_file.readline().split(' '))-1;
	feat_out.write(str(len(G.node)) + " " + str(feat_count) + "\n")
	count = 0

	feat_file = open(feat_filename)
	for lines in feat_file:
		words = lines.rstrip('\n').split(' ')
		words.pop(0)
		word_count = 0;
		for word in words:
			word_count+=1
			if (word=='1'):
				feat_out.write(str(word_count)+' ')
		feat_out.write('\n')
		count = count + 1;
		if (count == len(G.node)):
			break;


# print(G.nodes())
# print(G.edges())

if __name__ == "__main__": main()

