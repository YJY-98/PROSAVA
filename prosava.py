# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 09:56:24 2019

@author: dugujian
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 07:47:18 2019

@author: dugujian
"""
import platform
import argparse
from shutil import rmtree
import os
import matplotlib.pyplot as plt
#from openpyxl import Workbook
import networkx as nx
from ete3 import Tree,TreeStyle
from ete3 import SeqMotifFace
from random import randint

def _read_args():
    """Set parameters to execute scripts from the command line"""
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-E',
        metavar = 'E-value',
        help='1E-3 or 1E-10',
    )
    parser.add_argument(
        '-T',
        metavar = 'data type',
        help='sequence or proteome',
    )
    parser.add_argument(
        '-I',
        metavar = 'path of input file',
        help='input file path (FASTA file)',
    )
    
    # parser.add_argument(
    #     '-R',
    #     metavar='root path',
    #     help='root path',
    # )
    parser.add_argument(
        '-M',
        metavar='Pfam model file path',
        help='Pfam model file path',
    )
    
    
    
    return parser.parse_args()
def _hmm_press(dbPath,rootPath,platForm):
    if(os.path.exists(dbPath[:-3]+'h3f') and os.path.exists(dbPath[:-3]+'h3i') and os.path.exists(dbPath[:-3]+'h3m') and os.path.exists(dbPath[:-3]+'h3p')):
        return 0
    else:
        if(platForm=='Linux'):
            os.system(os.path.join(rootPath,'Linux/binaries/hmmpress ')+dbPath)
            return True
        elif(platForm=='Darwin'):    
            os.system(os.path.join(rootPath,'Darwin/binaries/hmmpress ')+dbPath)
            return True
        else:
            print('Unsupported Operating System (only for Linux or macOS)')
            return False
def _hmm_scan(eValue,dbPath,proPath,hmmPath,rootPath,platForm):
    """Recognition of the operating system and use hmmer software package to predict protein domains using Pfam hidden Markov models.
    
    Parameters
    ----------
    
    dbPath:str
        pathway for Pfam hidden Markov models
    proPath:str
        pathway for protein sequences/proteome
    rootPath:str
        pathway for prosava script
    ----------
    """
        
    if(platForm=='Linux'):#Report sequences and domains with E-values ​​less than designated threshold in the results, save the parsing table of the hit results for each domain to a file
        os.system(os.path.join(rootPath,'Linux/binaries/hmmsearch ')+'-E '+eValue+' --domE '+eValue+' --domtblout '+hmmPath+' '+dbPath+' '+proPath)
        return True
    elif(platForm=='Darwin'):    
        os.system(os.path.join(rootPath,'Darwin/binaries/hmmsearch ')+'-E '+eValue+' --domE '+eValue+' --domtblout '+hmmPath+' '+dbPath+' '+proPath)
        return True
    else:
        print('Unsupported Operating System (only for Linux or macOS)')
        return False
    
def _data_store(hmmPath,dataDict,eValue):
    m=open(hmmPath,'r')
    for line in m.readlines():
        if(line.startswith('#')):            
            continue
        else:
            parts=line.split()
            if(float(parts[12])>float(eValue)):
                continue
            if(parts[0] not in dataDict):
                dataDict[parts[0]]=[]
                dataDict[parts[0]].append([parts[3],[int(parts[17]),float(parts[12]),parts[4]]])
            else:
                dataDict[parts[0]].append([parts[3],[int(parts[17]),float(parts[12]),parts[4]]])
def _data_range(dataDict):
    for pro in dataDict:
        for i in range(len(dataDict[pro])):
            for j in range(i+1,len(dataDict[pro])):
                if(dataDict[pro][i][1][0]>dataDict[pro][j][1][0]):
                    dataDict[pro][i],dataDict[pro][j]=dataDict[pro][j],dataDict[pro][i]                 
def _data_analyse(dataDict,domainDict):
    for pro in dataDict:
        domainDict[pro]=[]
        for domain in dataDict[pro]:
            domainDict[pro].append(domain[0])
def _data_type(domainDict,typeList):
    kind=[]
    for pro in domainDict:
        if(domainDict[pro] not in kind):
            kind.append(domainDict[pro])
            typeList.append([domainDict[pro],1])
        else:
            typeList[kind.index(domainDict[pro])][1]+=1
    for i in range(len(typeList)):
        for j in range(i+1,len(typeList)):
            if(typeList[i][1]<typeList[j][1]):
                typeList[i],typeList[j]=typeList[j],typeList[i]
def _seq_analysis(dataDict,finalPath):
    m=open(os.path.join(finalPath,'seqAnalysis.csv'),'w')
    for seq in dataDict:
        m.write(seq)
        for domain in dataDict[seq]:
            m.write(','+domain[0]+','+str(domain[1][1]))
        m.write('\n')
            
        
def _domain_visual(finalPath,typeList):
    """Visualize the data in the typeList list, and finally generate a domain distribution pattern diagram
    
    Parameters
    ----------
    finalPath：str or unicode
        pathway for results
    typeList：list
        domains storage
    colorMap:dict
        domain color storage
    typeName：str
        domain names
    ----------
    """
    maxNum=0
    for ty in typeList:
        if(len(ty[0])>maxNum):
            maxNum=len(ty[0])
   
    plt.figure(figsize=(3*maxNum, len(typeList)*0.5))#size (length and width)
    colorMap={}
    for tl in typeList:
        for domain in tl[0]:
            if(domain not in colorMap):
                colorMap[domain]=''#Store all protein domains in the colorMap dictionary
    
    for domain in colorMap:
        colorMap[domain]=randomcolor()

    
    for tp in range(len(typeList)):
        plt.barh(left=0, y=len(typeList)-tp-1, width=len(typeList[tp][0])*120, color='0.4', height=0.3, orientation='horizontal')
        left=-110
        for index in range(len(typeList[tp][0])):
            left+=120
            plt.barh(left=left, y=len(typeList)-tp-1, width=100, color=colorMap[typeList[tp][0][index]], height=0.5, orientation='horizontal', alpha=1)
            plt.text(left,len(typeList)-tp-1.1,typeList[tp][0][index],color='w',alpha=1)
    typeName=['Type '+str(len(typeList)-i)+' ['+str(typeList[len(typeList)-i-1][1])+']' for i in range(len(typeList))]
    plt.yticks(range(len(typeList)),typeName, rotation=0,fontproperties = 'Times New Roman',size=12)
    plt.xticks([])
    plt.savefig(os.path.join(finalPath,'domainPic.pdf'), dpi=600)
def _domain_distribution(typeList,finalPath):
    
    """Output a table that stores domain information"""
    m=open(os.path.join(finalPath,'domainDistribution.csv'),'w')
    for ty in typeList:
        m.write(str(ty[1]))
        for domain in ty[0]:
            m.write(','+domain)
        m.write('\n')
    m.close()
    
def _domain_count(typeList,domainCount):
    for typ in typeList:
        for domain in typ[0]:
            if(domain not in domainCount):
                domainCount[domain]=typ[1]
            else:
                domainCount[domain]+=typ[1]
def _domain_frequence(domainCount,finalPath,domainFrequence):
    values=domainCount.values()
    index=list(set(values))
    index.sort()
    for domain in domainCount:
        if domainCount[domain] not in domainFrequence:
           domainFrequence[domainCount[domain]]=[]
           domainFrequence[domainCount[domain]].append(domain)
        else:
           domainFrequence[domainCount[domain]].append(domain)
    m=open(os.path.join(finalPath,'Domain_Frequence.csv'),'w')
    for i in index:
        m.write(str(i))
        for domain in domainFrequence[i]:
            m.write(','+domain)
        m.write('\n')
    m.close()
def _domain_interact(finalPath,typeList,domainCount,hubNodes,nodeInteract,nodeDegree):
    plt.figure(figsize=(50, 50))
    G = nx.Graph()

    allNum=0  #Total number of sequences
    for typ in typeList:
        allNum+=typ[1]

    for typ in typeList:
        if(len(typ[0])==1):
            continue
        else:
            domain=typ[0]
            for j in range(len(domain)-1):
                for k in range(j+1,len(domain)):
                    if(domain[j]==domain[k]):
                        continue
                    value=typ[1]
                    node_exist=False
                    for node in nodeInteract:
                        if(domain[j] in node and domain[k] in node):
                            node[2]+=value
                            node_exist=True
                    if(node_exist==False):
                        nodeInteract.append([domain[j],domain[k],value])
                        
#    for typ in typeList:
#        if(len(typ[0])==1):
#            continue
#        else:
#            domain=typ[0]
#            for j in range(len(domain)-1):
#                for k in range(j+1,len(domain)):
#                    if(domain[j]==domain[k]):
#                        continue
#                    value=round((1-0.2*(k-j-1))*(1+typ[1]/allNum),2)
#                    if(value<0):
#                        value=0
#                    node_exist=False
#                    for node in node_weight:
#                        if(domain[j] in node and domain[k] in node):
#                            node[2]=round(node[2]+value,2)
#                            node_exist=True
#                    if(node_exist==False):
#                        node_weight.append([domain[j],domain[k],value])

#    nodeInter=[i[2] for i in nodeInteract]
#    minInter=min(nodeInter)
#    maxInter=max(nodeInter)
#    addItem=9.5/(maxInter-minInter)
    G.add_weighted_edges_from(nodeInteract) 
    components = nx.connected_components(G)
#    print('num of connected_components:', nx.number_connected_components(G))
    largest_component = max(components, key=len)
    subgraph = G.subgraph(largest_component) #subgraph.edges and subgraph.nodes represent the number of nodes and edges in the subgraphs
    colors=[]
    node_size=[[],[]]
    
    for node in subgraph.nodes():
        if(subgraph.degree(node)>=30):
            colors.append('#FF6666')
            
            hubNodes.append(node)
            node_size[0].append(node)
            node_size[1].append(subgraph.degree(node))
        else:
            colors.append('#6699FF')
            node_size[0].append(node)
            node_size[1].append(subgraph.degree(node))
        
    for node in G.nodes():
        nodeDegree[node]=G.degree(node)
    #max_node=max(node_size[1])
    nodes_size=[i*50 for i in node_size[1]]
    subgraph_all=0
    for node in subgraph.nodes():
        subgraph_all+=domainCount[node]
    m=open(os.path.join(finalPath,'networkAppendix.txt'),'w')
    m.write('Total domain: '+str(sum(domainCount.values()))+'\n')
    m.write('All edges: '+str(sum([i[2] for i in nodeInteract]))+'\n')
    m.write('Non-redundant domain: '+str(len(G.nodes))+'\n')
    m.write('Maximum domain graph\n')
    m.write('No. of Nodes: '+str(len(subgraph.nodes))+'\n')
    m.write('No. of Edges: '+str(len(subgraph.edges))+'\n')
    m.write('Names of Hub Nodes(Occurrence>=30 times):')
    for node in hubNodes:
        m.write('\n'+node)
    m.close()
   
    pos=nx.spring_layout(subgraph)
    #nx.draw(subgraph,pos,with_labels=True,node_color=colors,node_size=nodes_size,font_size=6,width=[(int(d['weight'])-minInter)*addItem+0.5 for (u,v,d) in subgraph.edges(data=True)])
    nx.draw(subgraph,pos,with_labels=True,node_color=colors,node_size=nodes_size,font_size=6,width=[int(d['weight'])*0.5 for (u,v,d) in subgraph.edges(data=True)])
    #nx.draw_networkx_edges(subgraph,pos,width=[float(d['weight']*1) for (u,v,d) in subgraph.edges(data=True)])
    plt.savefig(os.path.join(finalPath,'domainInteract.jpg'), transparent=True,dpi=200)         
#    pos=nx.spring_layout(G)
#    plt.figure(figsize=(50, 50))
#    nx.draw(G,pos)
#    plt.savefig(finalPath+'/'+'domainInteract1.pdf', transparent=True,dpi=600)       
def _hub_table(hubNodes,pfDict,dataDict,finalPath,nodeDegree):
    for sq in dataDict:
        for domain in dataDict[sq]:
            if(domain[0] not in pfDict):
                pfDict[domain[0]]=domain[1][2]
    m=open(os.path.join(finalPath,'Hub_Nodes.csv'),'w')
    m.write('Domain,Pfam ID,Interactions\n')
    for node in hubNodes:
        m.write(node+','+pfDict[node]+','+str(nodeDegree[node])+'\n')
    m.close()
def _find_maxevalue(dataDict,node):
    maxEvalue=0
    for seq in dataDict:
        for domain in dataDict[seq]:
            if(domain[0]!=node):
                continue
            else:
                if(domain[1][1]>maxEvalue):
                    maxEvalue=domain[1][1]
                else:
                    continue
    return maxEvalue
def _find_minevalue(dataDict,node):  
    minEvalue=1
    for seq in dataDict:
        for domain in dataDict[seq]:
            if(domain[0]!=node):
                continue
            else:
                if(domain[1][1]<minEvalue):
                    minEvalue=domain[1][1]
                else:
                    continue
    return minEvalue
def _cyto_csv(finalPath,typeList,nodeInteract,nodeDegree,dataDict):
    m=open(os.path.join(finalPath,'Domain_Interactions.csv'),'w')
    
    for interact in nodeInteract:
        m.write(interact[0]+','+interact[1]+','+str(interact[2])+'\n')
    m.close()
    m=open(os.path.join(finalPath,'Non-redundant_Domains.csv'),'w')
    for node in nodeDegree:
        m.write(node+','+str(nodeDegree[node])+','+str(_find_maxevalue(dataDict,node))+','+str(_find_minevalue(dataDict,node))+'\n')
    m.close()
def _typeList_Filter(typeList,typeListFilter): 
    h = 10 if len(typeList)>10 else len(typeList)
    for tp in range(h):
        typeListFilter.append(typeList[tp])
def _pick_seq(typeListFilter,domainDict,finalPath,proPath,seqDict,dataDict_ty,cachePath):
    sequence=os.path.join(cachePath,'1_sequence')
    alignment=os.path.join(cachePath,'2_alignment')
    hmm=os.path.join(cachePath,'3_hmm_model')
    os.mkdir(sequence)
    os.mkdir(alignment)
    os.mkdir(hmm)
    tyList=[i[0] for i in typeListFilter]
    for i in range(len(tyList)):
        seqDict[i]=[]
    for seq in domainDict:
        try:
            dataDict_ty[seq]=tyList.index(domainDict[seq])
        except:
            continue
    for line in open(proPath,'r').readlines():
        if(line.startswith('>')):
            seqName=line.split()[0][1:]
            try:
                seqDict[dataDict_ty[seqName]].append(line)
            except:
                continue
        else:
            try:
                seqDict[dataDict_ty[seqName]][len(seqDict[dataDict_ty[seqName]])-1]+=line
            except:
                continue
    for tp_index in seqDict:
        m=open(os.path.join(sequence,'Type_'+str(tp_index)),'w')
        for seq in seqDict[tp_index]:
            m.write(seq)
        m.close()
def _seq_align(rootPath,cachePath,platForm):
    for file in os.listdir(os.path.join(cachePath,'1_sequence')):
        if(platForm=='Linux'):#在结果中报告E值小于1e-10的序列和domain，保存每个结构域的命中结果的解析表到文件中
            os.system(os.path.join(rootPath,'muscle_linux64 ')+'-in '+os.path.join(cachePath,'1_sequence',file)+' -out '+os.path.join(cachePath,'2_alignment',file)+' -maxiters 1 -diags -sv -distance1 kbit20_3')
        
        elif(platForm=='Darwin'):    
            os.system(os.path.join(rootPath,'muscle_darwin64 ')+'-in '+os.path.join(cachePath,'1_sequence',file)+' -out '+os.path.join(cachePath,'2_alignment',file)+' -maxiters 1 -diags -sv -distance1 kbit20_3')
       
        else:
            print('Unsupported Operating System (only for Linux or macOS)')
            exit()
def _hmm_build(rootPath,cachePath,platForm):
    for file in os.listdir(os.path.join(cachePath,'2_alignment')):
        if(platForm=='Linux'):#在结果中报告E值小于1e-10的序列和domain，保存每个结构域的命中结果的解析表到文件中
            os.system(os.path.join(rootPath,'Linux/binaries/hmmbuild ')+os.path.join(cachePath,'3_hmm_model',file+'.hmm ')+os.path.join(cachePath,'2_alignment',file))
        
        elif(platForm=='Darwin'):            
            os.system(os.path.join(rootPath,'Darwin/binaries/hmmbuild ')+os.path.join(cachePath,'3_hmm_model',file+'.hmm ')+os.path.join(cachePath,'2_alignment',file))
                
        else:
            print('Unsupported Operating System (only for Linux or macOS)')
            exit()
def _hmm_emit(cachePath,typeListFilter,rootPath,platForm):
    m=open(os.path.join(cachePath,'sequence.fasta'),'w')
    for file in os.listdir(os.path.join(cachePath,'3_hmm_model')):
        if(platForm=='Linux'):
            a=os.popen(os.path.join(rootPath,'Linux/binaries/hmmemit ')+os.path.join(cachePath,'3_hmm_model',file))
        elif(platForm=='Darwin'):
            a=os.popen(os.path.join(rootPath,'Darwin/binaries/hmmemit ')+os.path.join(cachePath,'3_hmm_model',file))
        else:
            print('Unsupported Operating System (only for Linux or macOS)')
            exit()
        for line in a.readlines():
            if(line.startswith('>')):
                index=int(file[5])
                #m.write('>Type'+str(index+1)+' （'+str(typeList[index][1])+'）\n')
                m.write('>Type '+str(index+1)+' ['+str(typeListFilter[index][1])+']\n')
            else:
                m.write(line)
    m.close()
def _muscle_handle(cachePath,platForm):
    if(platForm=='Linux'):
        os.system(os.path.join(rootPath,'muscle_linux64 ')+'-in '+os.path.join(cachePath,'sequence.fasta ')+'-out '+os.path.join(cachePath,'sequence.fasta.aln'))
        os.system(os.path.join(rootPath,'muscle_linux64 ')+'-maketree -in '+os.path.join(cachePath,'sequence.fasta.aln ')+'-out '+os.path.join(cachePath,'sequence.fasta.phy ')+'-cluster neighborjoining')
    elif(platForm=='Darwin'):
        os.system(os.path.join(rootPath,'muscle_darwin64 ')+'-in '+os.path.join(cachePath,'sequence.fasta ')+'-out '+os.path.join(cachePath,'sequence.fasta.aln'))
        os.system(os.path.join(rootPath,'muscle_darwin64 ')+'-maketree -in '+os.path.join(cachePath,'sequence.fasta.aln ')+'-out '+os.path.join(cachePath,'sequence.fasta.phy ')+'-cluster neighborjoining')
    else:
        print('Unsupported Operating System (only for Linux or macOS)')
        exit()
   
    if(os.path.exists(os.path.join(cachePath,'sequence.fasta.phy'))):
        return True
    else:
        return False
def _motif_datas(typeListFilter,motifData):
    colorMap={}
    for tl in typeListFilter: 
        for domain in tl[0]:
            if(domain not in colorMap):
                colorMap[domain]=''

    for domain in colorMap:
        colorMap[domain]=randomcolor()
        
    for i in range(len(typeListFilter)):
        motifData[i]=[]
        left=-50
        for j in typeListFilter[i][0]:

            left+=60
            motifData[i].append([left, left+50, "[]", None, 10, "black", colorMap[j], "arial|3|black|"+j])
def _tree_visual(typeListFilter,motifData,hmmPath,finalPath,cachePath):
    t = Tree(os.path.join(cachePath,'sequence.fasta.phy'))
    ts = TreeStyle()


    

    for sq in range(len(typeListFilter)):
        seqFace = SeqMotifFace((10+60*len(motifData[sq]))*'A', motifs=motifData[sq],seq_format="-")
        (t & 'Type '+str(sq+1)+' ['+str(typeList[sq][1])+']').add_face(seqFace, 0, "aligned")
    t.render(os.path.join(finalPath,'phylogeneticTree.pdf'), tree_style=ts)
def _judge_press(dbPath,platform,rootPath):

    if(os.path.isfile(dbPath[:-3]+'h3f')):
        print('Model file formatting Done')
    else:
        print('Formatting,please waiting for a moment!')
        if(platForm=='Linux'):#Report sequences and domains with E-values ​​less than designated threshold in the results, save the parsing table of the hit results for each domain to a file
            os.system(os.path.join(rootPath,'Linux/binaries/hmmpress ')+dbPath)
      
        elif(platForm=='Darwin'):    
            os.system(os.path.join(rootPath,'Darwin/binaries/hmmpress ')+dbPath)
      
        
        
def randomcolor():
    colorArr = ['9','A','B','C','D','E']
    color = ""
    for i in range(6):
        color += colorArr[randint(0,5)]
    return "#"+color
if __name__ == '__main__':
    
    args = _read_args()
    #rootPath=args.R
    rootPath = os.path.split(os.path.realpath(__file__))[0] #Absolute pathway
    try:
        os.mkdir(os.path.join(rootPath,'domain_dissection'))
    except:
        rmtree(os.path.join(rootPath,'domain_dissection'))
        os.mkdir(os.path.join(rootPath,'domain_dissection'))
    
    try:
        os.mkdir(os.path.join(rootPath,'result'))
    except:
        rmtree(os.path.join(rootPath,'result'))
        os.mkdir(os.path.join(rootPath,'result'))
    eValue=args.E  
    proPath=args.I  #pathway for sequence file
    finalPath=os.path.join(rootPath,'result')  #pathway for result output
    seqType=args.T #sequence type (protein seqeunces or proteome)
    dbPath=args.M  #pathway for pfam hmm models
    platForm=platform.system()
    _judge_press(dbPath,platForm,rootPath)
    
    cachePath=os.path.join(rootPath,'phylotree')  #pathway to the intermediate (temporary) files
    hmmPath =os.path.join(rootPath,'domain_dissection',os.path.splitext(os.path.split(proPath)[1])[0]+'.txt') #hmm result pathway

    
    
    dataDict={} #Sequence dictionary containing a list of domain-related parameters
    domainDict={} #Sequence dictionary with domain names for visualization
    typeList=[] #Cluster list, where elements contain domain information and the number of domain distribution patterns
    typeListFilter=[]
    typeSeq=[]
    seqDict={}
    motifData={}
    domainCount={} #Domain number dictionary, which stores the number of occurrences of each domain
    dataDict_ty={}
    domainFrequence={}
    hubNodes=[] #Storage of domains with equal to or more than 30 connectionns in the network
    pfDict={} #Pfam ID of domains in dictionary
    nodeDegree={} #
    nodeInteract=[] #Store protein domain association information, including protein domains and their association times
    
    _hmm_press(dbPath,rootPath,platForm)
    dealtResult=_hmm_scan(eValue,dbPath,proPath,hmmPath,rootPath,platForm)
    if(os.path.exists(hmmPath) and dealtResult):
        print('analysing please hang on......')
    else:
        exit()
    _data_store(hmmPath,dataDict,eValue) #Parse the domain data table and store it in the data dictionary
    _data_range(dataDict) #Sorting the positions of domains in proteins based on position information
    _data_analyse(dataDict,domainDict) #Simplified sequence dictionary with only domain name information
    _data_type(domainDict,typeList) #Clustering by domain architecture pattern
    
    if(seqType=='proteome'): #Perform the following analysis when the type of sequence is proteome
        try:
            rmtree(os.path.join(rootPath,'phylotree'))
        except:
            print('')
        _domain_count(typeList,domainCount)  #Count occurrences of each domain
        
        _domain_interact(finalPath,typeList,domainCount,hubNodes,nodeInteract,nodeDegree) #Visualize the domain-associated network and output related parameters to describe the network
        _cyto_csv(finalPath,typeList,nodeInteract,nodeDegree,dataDict) #Output table for CytoScape visualization
        _domain_frequence(domainCount,finalPath,domainFrequence)
        _hub_table(hubNodes,pfDict,dataDict,finalPath,nodeDegree)
    if(seqType=='sequence'):
        try:
            os.mkdir(os.path.join(rootPath,'phylotree'))
        except:
            rmtree(os.path.join(rootPath,'phylotree'))
            os.mkdir(os.path.join(rootPath,'phylotree'))
        _seq_analysis(dataDict,finalPath) #Output Uniprot ID and all domains and evalues
        
        _domain_visual(finalPath,typeList) #Distribution visualization
        _domain_distribution(typeList,finalPath) #Distribution table
        _typeList_Filter(typeList,typeListFilter) #Top 10 in typeList
        _pick_seq(typeListFilter,domainDict,finalPath,proPath,seqDict,dataDict_ty,cachePath)
        _seq_align(rootPath,cachePath,platForm)
        _hmm_build(rootPath,cachePath,platForm)
        
        _hmm_emit(cachePath,typeListFilter,rootPath,platForm)
        result=_muscle_handle(cachePath,platForm)
        if(result):
            print('analysing please hang on......')
        else:
            exit()
        _motif_datas(typeListFilter,motifData)
        _tree_visual(typeListFilter,motifData,hmmPath,finalPath,cachePath)
    
    #_domain_excel_output(typeList,finalPath)
    
    
    
    
   
    #words = ' '.join(args.words)
#    if args.db is not None:
#        words = '%s: %s' % (args.db, args.words)
   # print(dbPath,proPath,finalPath)
