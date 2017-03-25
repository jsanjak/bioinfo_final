import sys
import os
import fnmatch

def makedir(base,sample,experiment):
    home=os.getcwd()
    os.chdir(base)
    if (not os.path.exists(sample)) :
        os.mkdir(sample)
    os.chdir(sample)
    if (not os.path.exists(experiment)):
        os.mkdir(experiment)
    os.chdir(home)
    path = "/".join([base,sample,experiment])
    return(path)

def main():
    home=os.getcwd()
    #DNA seq stuff
    
    #RNA seq stuff
    path = "rawdata"
    finfo = open(path+"/RNAseq384_SampleCoding.txt","r")
    RNAseq_info = finfo.readlines()
    finfo.close()
    
    tissues = {"BodyPlate":"Body","HeadPlate":"Head","EmbyroPlate":"Embryo","PupaPlate":"Pupa"}
    home=os.getcwd()
    for i in range(len(RNAseq_info)-1):
        for flowcell in range(1,5,1):
            inf = RNAseq_info[i+1].split()
            sample = inf[0]
            plex = inf[1].split("_")[0]
            tiss = tissues[inf[4]]
            samplepath = path + "/RNAseq384plex_flowcell0"+str(flowcell)+"/Project_" + plex+ "/Sample_" + inf[0] +"/"  
            RNA_seqs = sorted(fnmatch.filter(os.listdir(samplepath),"*.fastq.gz"))
            #make the directory for the symlinked data
            os.chdir("data/RNAseq")
            if (not os.path.exists(tiss)):
                os.mkdir(tiss)           
            os.chdir(home)
            fpath = "/".join(["data/RNAseq",tiss])
            fnameF = fpath + "/Flowcell0"+str(flowcell)+".Sample"+inf[0]+".RIL." + inf[8] + "."+ inf[3] +  ".F.fq.gz"
            fnameR = fpath + "/Flowcell0"+str(flowcell)+".Sample"+inf[0]+".RIL." + inf[8] + "."+ inf[3] +  ".R.fq.gz"
            if ( not os.path.islink(fnameR) ):
                os.symlink(home + "/" + samplepath + RNA_seqs[1], fnameR)
            if ( not os.path.islink(fnameF) ):
                os.symlink(home + "/" + samplepath + RNA_seqs[0], fnameF)
     
if __name__ == "__main__":
    main()
