import networkx as nx
import numpy as np
import itertools
from numpy import linalg as LA
import math
from operator import itemgetter

def loadGraph(edge, feat):
	#Store feat vec as ndarray
	featMap = {}
	lines = []
	with open(feat) as f:
		lines = f.readlines()
		lines = map(lambda x : map(int, x.split()), lines)
	for line in lines:
		featMap[line[0]] = np.ndarray(shape=(len(line)-1, 1), dtype=int, buffer=np.array(line[1:]))
	
	G = nx.read_edgelist(edge, nodetype=int)
	return G, featMap

# the first element of each line of community is by default the name of the community
# community is a list of list containing int id
def loadCommunity(communityfile):
	lines = []
	with open(communityfile) as f:
		lines = f.readlines()
	lines = map(lambda x : x.split(), lines)
	for line in lines:
		line[1:] = map(int, line[1:])
	return lines

def calculateDensity(Graph, community):
	result = []
	for com in community:
		subg = Graph.subgraph(com[1:])
		# print subg.nodes()
		result.append(nx.density(subg))
	return result

def calculateHomo_meanNormalized(Graph, nodefeat, community):
	result = []
	for com in community:
		subg = Graph.subgraph(com[1:])
		nodes = subg.nodes()
		sz = len(nodes)
		homo = 0.0
		var = np.var([LA.norm((nodefeat[x[0]] -  nodefeat[x[1]]), 2) for x in list(itertools.combinations(nodes, 2))])
		if sz > 1:
			for pair in itertools.combinations(nodes, 2):
				homo += math.exp((-1 * LA.norm((nodefeat[pair[0]] - nodefeat[pair[1]]), 2))/var)
			homo /= (sz * (sz - 1) / 2)
		result.append(homo)
	return result


# This is calculated using the linear interpolation of most important features in the community
# The metric is calculated using weights times the transpose of important features of same size
# Sum of weights has to equal 1.0 to make sense
def calculateHomo_linearInterpolation(Graph, nodefeat, community, weights):
	num_feat = len(weights)
	result = []
	for com in community:
		subg = Graph.subgraph(com[1:])
		nodes = subg.nodes()
		if len(nodes) > 0:
			cumulated_feat = reduce(lambda x,y : x+y, [nodefeat[x] for x in nodes if nodefeat[x] is not None]).tolist()
			cumulated_feat = sorted(zip(cumulated_feat, range(len(cumulated_feat))), key=itemgetter(0), reverse=True)
			top_k_feat = [x[0][0]/float(len(nodes)) for x in cumulated_feat[:len(weights)]]
			homo = np.dot(np.asarray(weights).reshape(1, len(weights)), np.asarray(top_k_feat).reshape(len(weights), 1))
			homo = homo[0][0]
		else:
			homo = 0
		result.append(homo)
	return result


def pairwiseDense(Graph, community):
	result = []
	for com in community:
		subg = Graph.subgraph(com[1:])
		subgsz = len(subg.nodes())
		try:
			tmp = [subg.degree(node)/float(subgsz-1) if subgsz > 1 else 0.0 for node in com[1:]]
		except nx.exception.NetworkXError:
			# this error means some nodes in the community is not in the subgraph
			tmp = [0.0]
		result.append(tmp)
	return result

def nodeHomoDiff(vec1, vec2, variance=1.0):
	return math.exp((-1 * LA.norm((vec1 - vec2), 2))/variance)


def pairwiseHomo_meanNormalized(Graph, nodefeat, community):
	result = []
	for com in community:
		subg = Graph.subgraph(com[1:])
		nodes = subg.nodes()
		sz = len(nodes)
		tmp = []
		var = np.var([LA.norm((nodefeat[x[0]] -  nodefeat[x[1]]), 2) for x in list(itertools.combinations(nodes, 2))])
		for node in com[1:]:
			try:
				cur = sum([nodeHomoDiff(nodefeat[node], nodefeat[other], variance=var) for other in subg.neighbors(node)])
				cur /= float(sz-1)
			except:
				cur = 0.0
			tmp.append(cur)
		result.append(tmp)
	return result


def plot_preliminary(folder, raw_edge, raw_feature, community_output):
	import matplotlib.pyplot as plt

	density_lis, homo_lis = [], []

	for (edgepath, featpath, case) in zip(raw_edge, raw_feature, community_output):
		G, featMap = loadGraph(edgepath, featpath)
		for compath in case:
			com = loadCommunity(compath)
			print compath
			fileName = compath.split("/")[-1]
			
			lis = calculateDensity(G,com)	
			lis = filter(lambda x : x > 0, lis)
			density_lis.append(lis)
			plt.hist(lis, 20, alpha=0.5)
			plt.xlabel('Desnity of the Community')
			plt.xlim([0,1])
			plt.ylabel('Number of Communities')
			plt.title(r'Histogram of Community density using %s' % fileName)
			plt.savefig(folder + '/%s_density_hist.png' % fileName)
			plt.clf()

			lis = calculateHomo_meanNormalized(G, featMap, com)
			lis = filter(lambda x : not math.isnan(x), lis)
			homo_lis.append(lis)
			plt.hist(lis, 20, alpha=0.5)
			plt.xlabel('Homogeneity of the Community')
			plt.ylabel('Number of Communities')
			plt.title(r'Histogram of Community Homogeneity using %s' % fileName)
			plt.savefig(folder + '/%s_homogeneous_hist.png' % fileName)
			plt.clf()

			lis = pairwiseDense(G, com)
			plt.boxplot(lis, patch_artist=True)
			plt.xlabel('Communities ID')
			plt.ylabel('Point-wise Desnity')
			plt.title(r'Boxplot of Point-wise Community Density using %s' % fileName)
			plt.savefig(folder + '/%s_point_density_boxplot.png' % fileName)
			plt.clf()

			lis = pairwiseHomo_meanNormalized(G, featMap, com)
			plt.boxplot(lis, patch_artist=True)
			plt.xlabel('Communities ID')
			plt.ylabel('Point-wise Homogeneity')
			plt.title(r'Boxplot of Point-wise Community Homogeneity using %s' % fileName)
			plt.savefig(folder + '/%s_point_homogeneous_boxplot.png' % fileName)
			plt.clf()

			# print calculateHomo_linearInterpolation(G, featMap, com, [0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1])

			# print calculateDensity(G,com)
			# print calculateHomo(G, featMap, com)
			# print pairwiseDense(G, com)
			# print pairwiseHomo(G, featMap, com)
			print "--------------------------"
	return density_lis, homo_lis



def get_pairwise(folder, raw_edge, raw_feature, community_output):

	density_lis, homo_lis = [], []

	for (edgepath, featpath, case) in zip(raw_edge, raw_feature, community_output):
		G, featMap = loadGraph(edgepath, featpath)
		com = loadCommunity(case)
		print case
		fileName = case.split("/")[-1]
		
		lis = pairwiseDense(G, com)
		lis = map(lambda x : map(str, x), lis)
		with open(folder + '/%s_point_density.txt' % fileName, 'w+') as f:
			f.write("\n".join([",".join(x) for x in lis]))
		# plt.boxplot(lis, patch_artist=True)
		# plt.xlabel('Communities ID')
		# plt.ylabel('Point-wise Desnity')
		# plt.title(r'Boxplot of Point-wise Community Density using %s' % fileName)
		# plt.savefig(folder + '/%s_point_density_boxplot.png' % fileName)
		# plt.clf()

		lis = pairwiseHomo_meanNormalized(G, featMap, com)
		lis = map(lambda x : map(str, x), lis)
		with open(folder + '/%s_point_homo.txt' % fileName, 'w+') as f:
			f.write("\n".join([",".join(x) for x in lis]))
		# plt.boxplot(lis, patch_artist=True)
		# plt.xlabel('Communities ID')
		# plt.ylabel('Point-wise Homogeneity')
		# plt.title(r'Boxplot of Point-wise Community Homogeneity using %s' % fileName)
		# plt.savefig(folder + '/%s_point_homogeneous_boxplot.png' % fileName)
		# plt.clf()

		# print calculateHomo_linearInterpolation(G, featMap, com, [0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1])

		# print calculateDensity(G,com)
		# print calculateHomo(G, featMap, com)
		# print pairwiseDense(G, com)
		# print pairwiseHomo(G, featMap, com)
		print "--------------------------"
	

def jaccard_sim(lis1, lis2):
	numerator = filter(lambda x : x in lis2, lis1)
	denumerator = list(set(lis1).union(lis2))
	return len(numerator)/float(len(denumerator))

def f1_sim(truth, com):
	try:
		correct_nodes = filter(lambda x : x in truth, com)
		precision = len(correct_nodes)/float(len(com))
		recall = len(correct_nodes)/float(len(truth))
		f1 = 2 * (precision * recall)/(precision + recall)
		return f1
	except:
		return 0
	# correct_nodes = filter(lambda x : x in truth, com)
	# precision = len(correct_nodes)/float(len(com))
	# recall = len(correct_nodes)/float(len(truth))
	# f1 = 2 * (precision * recall)/(precision + recall)

	return 0


def experiment_metric(truth, community):
	t_lis = [loadCommunity(x) for x in truth]
	c_lis = [loadCommunity(x) for x in community]
	f1_result, jaccard_result = [], []
	for (t, c) in zip(t_lis, c_lis):
		t = map(lambda x : x[1:], t)
		c = map(lambda x : x[1:], c)
		# do detected
		temp_result_f1_a, temp_result_jac_a = 0.0, 0.0
		for current_detected_com in c:
			temp_result_f1_a += max([f1_sim(current_detected_com, x) for x in t])
			temp_result_jac_a += max([jaccard_sim(current_detected_com, x) for x in t])

		# do truth
		temp_result_f1_b, temp_result_jac_b = 0.0, 0.0
		for current_truth_com in t:
			temp_result_f1_b += max([f1_sim(current_truth_com, x) for x in c])
			temp_result_jac_b += max([jaccard_sim(current_truth_com, x) for x in c])

		f1_value = temp_result_f1_a/(2*float(len(c))) + temp_result_f1_b/(2*float(len(t)))
		jac_value = temp_result_jac_a/(2*float(len(c))) + temp_result_jac_b/(2*float(len(t)))
		
		f1_result.append(f1_value)
		jaccard_result.append(jac_value)
	return f1_result, jaccard_result
	


def plot_experiment_experiment(truth, list_com, list_name):
	import matplotlib.pyplot as plt
	color_sequence = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c','#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5','#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f','#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']

	f1_list, jacc_list = [], []
	for rank, name in enumerate(list_name):
		print "processing %s" % (name)
		f1,jacc = experiment_metric(truth, list_com[rank])
		f1_list.append(f1)
		jacc_list.append(jacc)

	return f1_list, jacc_list
	


fb_circle_lis = map(lambda x : str(x),[0, 107, 348, 414, 686, 698, 1684, 1912, 3437, 3980])

fb_raw_edge = ["data/facebook/raw/%s.edges" % (num) for num in fb_circle_lis]
fb_raw_feature = ["data/facebook/raw/%s.feat" % (num) for num in fb_circle_lis]

# community_output = ["testdata/1912_truth.txt","testdata/1912_cesna.txt", "testdata/1912_bigcalm.txt", "testdata/1912_CNM.txt"]
truth_community_output = ["testdata/fb/%s_truth.txt" % (num) for num in fb_circle_lis]
cesna_community_output = ["testdata/fb/%s_cesna.txt" % (num) for num in fb_circle_lis]
bigclam_community_output = ["testdata/fb/%s_bigclam.txt" % (num) for num in fb_circle_lis]
circles_community_output = ["testdata/fb/%s_circles.txt" % (num) for num in fb_circle_lis]
mincut_community_output = ["testdata/fb/%s_mincut.txt" % (num) for num in fb_circle_lis]
infomap_community_output = ["testdata/fb/%s_infomap.txt" % (num) for num in fb_circle_lis]
community_output = map(lambda x : list(x), zip(truth_community_output, cesna_community_output, bigclam_community_output, circles_community_output, mincut_community_output, infomap_community_output))


'''
Uncomment this sections for plot all png into plots folder


a,b = plot_preliminary("plots", fb_raw_edge, fb_raw_feature, community_output)
print a
print b 

'''


'''
Uncomment this section for store all pointwise_data into pointwise_data folder


get_pairwise("pointwise_data", fb_raw_edge, fb_raw_feature, truth_community_output)
'''


'''
Uncomment this section for get evaluation performance, change the second paramter to other community ouput defined above if necessary


f1,jacc = experiment_metric(truth_community_output, cesna_community_output)
'''