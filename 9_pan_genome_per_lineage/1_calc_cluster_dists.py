import argparse
import os
from numpy import mean, median
import time

def assemblies_to_gffs(metadata_file):
    ''' parse the input metadata file so that each assembly
    has its equivalent GFF file.
    The GFF files are required for running ROARY.
    return: dictionary of assembly -> GFF'''
    print("Getting GFF filenames...")
    fa_to_gff = {}
    with open(metadata_file) as f:
        for line in f:
            toks = line.strip().split("\t")
            if toks[0] == "ID":
                assembly_index = toks.index("Assembly_Location")
                annot_index = toks.index("Annotation_Location")
                continue
            assemblies = toks[assembly_index].split(",")
            gffs = toks[annot_index].split(",")
            if len(gffs) == 1:
                for a in assemblies:
                    fa_to_gff[a] = gffs[0]
            else:
                for a in assemblies:
                    basename = a.split("/")[-1]
                    basename = basename.split(".")[0]
                    for g in gffs:
                        if basename in g:
                            fa_to_gff[a] = g
                    if a not in fa_to_gff: ## SANITY CHECK - SHOULDNT HAPPEN
                        print("NOT FOUND! \n assembly: %s\n gff: %s" % (str(a), str(gffs)))
    return fa_to_gff


def get_cluster_sizes(clusters_file, out):
    ''' go over the poppunk cluster output
    For each assembly, assign it a cluster membership
    Return 1) dictionary assembly -> cluster
        2) For each cluster, calculate its size. If it's smaller
    than min_cluster_size disregard it for ROARY by
    saving it in the cluster_size dictionary
    Output: a file with each clusters size '''
    print("Calculating cluster sizes...")
    try:
        os.makedirs(out)
    except Exception:
        pass
    clusters = {}
    cluster_size = {}
    with open(clusters_file) as f:
        for line in f:
            if line.startswith("Taxon"):
                continue
            toks = line.strip().split(",")
            if toks[1] not in cluster_size:
                cluster_size[toks[1]] = 0
            cluster_size[toks[1]] += 1
            clusters[toks[0]] = toks[1]
    ## generate cluster size output
    with open(os.path.join(out, "cluster_sizes.csv"),"w") as out:
        out.write("Cluster, Size\n")
        for cluster in cluster_size:
            out.write(cluster + "," + str(cluster_size[cluster]) + "\n")
    return (clusters, cluster_size)


def get_updated_cluster_sizes():
    cluster_sizes = {}
    with open("cluster_sizes_updated.csv") as f:
        for line in f:
            if line.startswith("Cluster"):
                continue
            toks = line.strip().split(",")
            cluster_sizes[toks[0]] = int(toks[1])
    return cluster_sizes

def metadata_per_cluster(metadata_file, cluster_sizes, clusters, min_cluster_size, out):
    ''' generate a new clusters output so that I can generate bar/piecharts
    where I know exactly how many of each ST, pathotype, etc. are in each of the
    poppunk clusters
    output: melted dataframe with cluster and count of each number from each category
    return: a dictionary with each of the metadata objects for each strains
    (to be used for the distances ouput)'''

    ## variables to ouput in file
    variables = ["Pathotype", "Publication", "Continent", "Year","Isolation", "ST", "MASH"]
    clusters_metadata = {}
    all_metadata = {}
    with open(metadata_file) as f:
        for line in f:
            toks = line.strip().split("\t")
            if line.startswith("ID"):
                indexs = [toks.index(i) for i in variables]
                assembly_index = toks.index("Assembly_Location")
                poppunk_index = toks.index("Poppunk_cluster")
                continue
            cluster = toks[poppunk_index]
            if toks[poppunk_index] not in cluster_sizes.keys():
                continue

            a = toks[assembly_index].split(",")[0]
            all_metadata[a] = {}
            if cluster not in clusters_metadata:
                clusters_metadata[cluster] = {}
                for v in variables:
                    clusters_metadata[cluster][v] = {}
            for i in range(len(variables)):
                curr_value = toks[indexs[i]]
                v = variables[i]
                ## fix the metadata outputs:
                curr_value = curr_value.replace("~","") ## remove "~" from STs
                if v == "Isolation" and curr_value not in ["feces","urine","blood"]:
                    curr_value = "other/unknown"
                elif "ehec" in curr_value:
                    curr_value = "ehec"
                elif "aeec" in curr_value:
                    curr_value = "epec/eaec"
                elif curr_value == "upec":
                    curr_value = "expec"
                elif v == "Year":
                    curr_value = curr_value.replace("s","")
                    curr_value = curr_value.split("-")
                    try:
                        curr_value = map(int, curr_value)
                        curr_value = mean(curr_value)
                    except Exception:
                        curr_value = "nd"
                elif "public health england" in curr_value:
                    curr_value = "phe"
                elif "fda" in curr_value or "food and drug administration" in curr_value:
                    curr_value = "fda"
                elif "sanger" in curr_value:
                    curr_value = "sanger"
                elif "centers for disease control" in curr_value or "cdc" in curr_value:
                    curr_value = "cdc"
                elif "kallonen2017_bsac" in curr_value:
                    curr_value = "bsac"
                elif "hazen2013" in curr_value:
                    curr_value = "hazen2013"
                elif "ingle2016" in curr_value:
                    curr_value = "ingle2016"

                ## this value has never been seen for this cluster
                if curr_value not in clusters_metadata[cluster][v]:
                    clusters_metadata[cluster][v][curr_value] = 0
                ## add count by one for this value in this variable in this cluster
                clusters_metadata[cluster][v][curr_value] += 1
                all_metadata[a][v] = curr_value

    ## generate the melted output file -> easy to work with in R
    out = open(os.path.join(out,"metadata_per_cluster.csv"),"w")
    out.write("cluster\tvariable\tvalue\tcount\n")
    for cluster in clusters_metadata:
        for var in clusters_metadata[cluster]:
            for value in clusters_metadata[cluster][var]:
                out.write("\t".join([cluster, var, str(value), str(clusters_metadata[cluster][var][value] / float(cluster_sizes[cluster]))]) + "\n")
    out.close()
    return all_metadata

def get_dist_within_cluster(clusters, cluster_sizes, dists_file, out, min_cluster_size, fa_to_gff, metadata):
    ''' Due to memory constrains, read the dists file line by line.
    Check membership of both clusters, if they are different ignore,
    if they are the same but the cluster is too small also ignore
    (i.e. skip MOST lines of any calculation)
    Otherwise, save the core and acc dists of this cluster
    to be used for ROARY input and save all the GFF files of this clusters to
    a file to be used by ROARY
    Return: dictionary with cluster -> core and accessory distance
    Output: a file summarising all the within-cluster distances with metadata'''
    print("Calculating within-cluster distances...")

    md_out = {}
    metadata_vals = metadata.values()[0].keys()
    header = ['{}_{}'.format(a, b) for b in ["1","2"] for a in metadata_vals]
    for cluster in cluster_sizes:
        if cluster_sizes[cluster] >= min_cluster_size:
            md_out[cluster] =  open(os.path.join(out, str(cluster) + "_metadata_within_cluster.csv"), "w")
            md_out[cluster].write("Core_dist, Acc_dist, " + ",".join(header) + "\n")
    within_cluster_dist = {}
    gffs = {}
    cnt = 0
    with open(dists_file) as f:
        for line in f:
            toks = line.strip().split("\t")
            if line.startswith("Query"):
                continue
            if clusters[toks[0]] != clusters[toks[1]] or cluster_sizes[clusters[toks[0]]] < min_cluster_size:
                continue
            cluster = clusters[toks[0]]
            if cluster not in within_cluster_dist:
                within_cluster_dist[cluster] = {"core":[], "acc":[]}
                gffs[cluster] = set()
            within_cluster_dist[cluster]["core"].append(float(toks[2]))
            within_cluster_dist[cluster]["acc"].append(float(toks[3]))
            gffs[cluster].add(fa_to_gff[toks[0]])
            gffs[cluster].add(fa_to_gff[toks[1]])

            md_out[cluster].write(",".join([toks[2], toks[3]]))
            for j in [0,1]: ## get all the metadata values for toks[0] and toks[1]
                for v in metadata_vals:
                    md_out[cluster].write("," + str(metadata[toks[j]][v]))
            md_out[cluster].write("\n")

            # if cnt == 10000: ## if testing uncomment here
            #     break
            # cnt += 1
    for cluster in md_out:
        md_out[cluster].close()
    ## calc the averages
    ## generate output so that this doesn't need to be repeated
    with open(os.path.join(out, "within_cluster_dist.csv"),"w") as f_out:
        f_out.write("Cluster, Core, Core_max, Core_median, Acc, Acc_max, Acc_median\n")
        for cluster in within_cluster_dist:
            core_max = max(within_cluster_dist[cluster]["core"])
            acc_max = max(within_cluster_dist[cluster]["acc"])
            mean_core = mean(within_cluster_dist[cluster]["core"])
            mean_acc =  mean(within_cluster_dist[cluster]["acc"])
            core_median = median(within_cluster_dist[cluster]["core"])
            acc_median = median(within_cluster_dist[cluster]["acc"])
            f_out.write(",".join(map(str, [cluster, mean_core, core_max, core_median,
                                            mean_acc, acc_max, acc_median])) + "\n")

    ## create the GFF outputs -> for running Roary
    out = os.path.abspath(os.path.join(out,"gff_jobs"))
    try:
        os.makedirs(out)
    except Exception:
        pass
    for cluster in gffs:
        with open(os.path.join(out, "jobs_" + cluster + ".txt"), "w") as f_out:
            for gff in gffs[cluster]:
                f_out.write(gff + "\n")
    return


def add_clusters_to_dists(between_cluster_dist, cluster1, cluster2, core, accessory):
    ''' add the distance between two clusters to the between clusters dictionary '''
    if cluster1 not in between_cluster_dist:
        between_cluster_dist[cluster1] = {}
    if cluster2 not in between_cluster_dist[cluster1]:
        between_cluster_dist[cluster1][cluster2] = {"core":[], "accessory":[]}
    between_cluster_dist[cluster1][cluster2]["core"].append(core)
    between_cluster_dist[cluster1][cluster2]["accessory"].append(accessory)
    return

def between_cluster_dist(clusters, cluster_sizes, min_cluster_size, dists_file, fa_to_gff, out, metadata):
    ''' read the dists file again, this time I only care about the
    distance between two clusters that I decided to keep, therefore
    the number of calculations on the whole file is still small '''
    print("Calculating between-cluster distances...")
    between_cluster_dist = {}

    md_out = open(os.path.join(out, "metadata_between_clusters.csv"), "w")
    metadata_vals = metadata.values()[0].keys()
    header = ['{}_{}'.format(a, b) for b in ["1","2"] for a in metadata_vals]
    md_out.write("Member_1, Member_2, Cluster_1, Cluster_2, Core_dist, Acc_dist, " + ",".join(header) + "\n")

    cnt = 0
    with open(dists_file) as f:
        for line in f:
            toks = line.strip().split("\t")
            if line.startswith("Query"):
                continue
            if  clusters[toks[0]] == clusters[toks[1]]: ## same cluster
                continue
            ## the cluster size is too small
            if cluster_sizes[clusters[toks[0]]] < min_cluster_size or cluster_sizes[clusters[toks[1]]] < min_cluster_size:
                continue
            cluster1 = int(clusters[toks[0]])
            cluster2 = int(clusters[toks[1]])
            if cluster1 < cluster2:
                add_clusters_to_dists(between_cluster_dist, cluster1, cluster2, float(toks[2]), float(toks[3]))
            else:
                add_clusters_to_dists(between_cluster_dist, cluster2, cluster1, float(toks[2]), float(toks[3]))
            md_out.write(",".join(map(str,[toks[0], toks[1], cluster1, cluster2, toks[2], toks[3]])))
            for j in [0,1]: ## get all the metadata values for toks[0] and toks[1]
                for v in metadata_vals:
                    md_out.write("," + str(metadata[toks[j]][v]))
            md_out.write("\n")
            # if cnt == 10000: ## if testing uncomment here
            #     break
            # cnt += 1
    md_out.close()

    cluster_ids = between_cluster_dist.keys()
    with open(os.path.join(out, "between_cluster_dist.csv"),"w") as f_out:
        f_out.write("cluster1, cluster2, core, core_max, core_median, accessory, accessory_max, accessory_median\n")
        for c1 in cluster_ids:
            for c2 in cluster_ids:
                if c2 <= c1:
                    continue
                core_max = max(between_cluster_dist[c1][c2]["core"])
                core_mean = mean(between_cluster_dist[c1][c2]["core"])
                core_median = median(between_cluster_dist[c1][c2]["core"])
                acc_max = max(between_cluster_dist[c1][c2]["accessory"])
                acc_mean = mean(between_cluster_dist[c1][c2]["accessory"])
                acc_median = median(between_cluster_dist[c1][c2]["accessory"])
                f_out.write(",".join(map(str, [c1, c2, core_mean, core_max, core_median, acc_mean, acc_max, acc_median])) + "\n")
    return


def run(args):
    args.out = os.path.abspath(args.out)
    args.metadata_file = os.path.abspath(args.metadata_file)
    fa_to_gff = assemblies_to_gffs(args.metadata_file)
    clusters, cluster_sizes = get_cluster_sizes(args.clusters_file, args.out)
    cluster_sizes = get_updated_cluster_sizes( )
    all_metadata = metadata_per_cluster(args.metadata_file, cluster_sizes, clusters, args.min_cluster_size, args.out)
    quit()
    within_cluster_dist = get_dist_within_cluster(clusters, cluster_sizes,
    args.dists_file, args.out, args.min_cluster_size, fa_to_gff, all_metadata)
    quit()
    between_cluster_dist(clusters, cluster_sizes, args.min_cluster_size,  args.dists_file, fa_to_gff, args.out, all_metadata)
    return


def get_options():
    parser = argparse.ArgumentParser(description='Calculate the between and within cluster distances. Create files with roary jobs to run.')
    # input options
    parser.add_argument('--clusters_file', required=True,
                        type=str,
                        help='poppunk clusters output file')
    parser.add_argument('--dists_file',
                        required=True,
                        type=str,
                        help='output of extract_distances.py of core and accessory distances between every two strains')
    parser.add_argument('--metadata_file',
                        required=True,
                        type=str,
                        help='Metadata file containing the GFF files of the assemblies')
    parser.add_argument('--out', help='Name of output directory [%(default)s]', default = "dists_analysis",
                        type=str,
                        required=False)
    parser.add_argument('--min_cluster_size', help='Minimum size of a cluster to be used in further analysis [%(default)s]',
                        type=int,
                        default=30)
   	# parser.add_argument('--merge_files', action='store_true', help='Merge the distances and clusters file into one', default=False)
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    start = time.time()
    options = get_options()
    run(options)
    end = time.time()
    print("Time: %s" %str(end - start))
