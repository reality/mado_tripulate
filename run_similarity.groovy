@Grab(group='org.semanticweb.elk', module='elk-owlapi', version='0.4.3')
@Grab(group='net.sourceforge.owlapi', module='owlapi-api', version='4.2.5')
@Grab(group='net.sourceforge.owlapi', module='owlapi-apibinding', version='4.2.5')
@Grab(group='net.sourceforge.owlapi', module='owlapi-impl', version='4.2.5')
@Grab(group='com.github.sharispe', module='slib-sml', version='0.9.1')

import groovyx.gpars.*
import org.codehaus.gpars.*
import java.util.concurrent.*

import org.openrdf.model.URI;
import slib.graph.algo.accessor.GraphAccessor;
import slib.graph.algo.extraction.utils.GAction;
import slib.graph.algo.extraction.utils.GActionType;
import slib.graph.algo.validator.dag.ValidatorDAG;
import slib.graph.io.conf.GDataConf;
import slib.graph.io.conf.GraphConf;
import slib.graph.io.loader.GraphLoaderGeneric;
import slib.graph.io.util.GFormat;
import slib.graph.model.graph.G;
import slib.graph.model.impl.graph.memory.GraphMemory;
import slib.graph.model.impl.repo.URIFactoryMemory;
import slib.graph.model.repo.URIFactory;
import slib.sml.sm.core.engine.SM_Engine;
import slib.sml.sm.core.metrics.ic.utils.IC_Conf_Topo;
import slib.sml.sm.core.metrics.ic.utils.ICconf;
import slib.graph.algo.extraction.utils.*
import slib.sglib.io.loader.*
import slib.sml.sm.core.metrics.ic.utils.*
import slib.sml.sm.core.utils.*
import slib.sglib.io.loader.bio.obo.*
import org.openrdf.model.URI
import slib.graph.algo.extraction.rvf.instances.*
import slib.sglib.algo.graph.utils.*
import slib.utils.impl.Timer
import slib.graph.algo.extraction.utils.*
import slib.graph.model.graph.*
import slib.graph.model.repo.*
import slib.graph.model.impl.graph.memory.*
import slib.sml.sm.core.engine.*
import slib.graph.io.conf.*
import slib.graph.model.impl.graph.elements.*
import slib.graph.algo.extraction.rvf.instances.impl.*
import slib.graph.model.impl.repo.*
import slib.graph.io.util.*
import slib.graph.io.loader.*

import slib.sml.sm.core.metrics.ic.utils.IC_Conf_Corpus;
import slib.sml.sm.core.utils.SMConstants;
import slib.sml.sm.core.utils.SMconf;
import slib.utils.ex.SLIB_Exception;
import slib.utils.impl.Timer;

def aFile = args[0]
def rmneg = args[1] == 'true'
def additionalAssocation = args[2] 
def oFile = args[3]

def patient_visit_diagnoses = [:]


//new File('./reconst_train.csv').splitEachLine('\t') {
//new File('./recon1k.csv').splitEachLine('\t') {
//new File('./reconst_train.csv').splitEachLine('\t') {
new File('./TEST_sampled_patient_visits.csv').splitEachLine('\t') {
  patient_visit_diagnoses[it[0]] = it[1].tokenize(',').collect { o -> 
  
    o = o.replaceAll('"','')
    o
  }
}

ConcurrentHashMap aMap = [:]

def cList = []

def aFileContent = []
new File(aFile).splitEachLine('\t') { aFileContent << it }

def t = 0
aFileContent.each {
  if(it[0] && it[1]) {
  it[0] = it[0].tokenize('.')[0]
  if(patient_visit_diagnoses.containsKey(it[0])) {
    if(rmneg && it[3] && it[3].indexOf('neg') != -1 ) { return ;}
    if(!aMap.containsKey(it[0])) {
      aMap[it[0]] = []
    }
    it[1] = it[1].replace('<', '').replace('>', '').tokenize('/').last()
    def z = it[1].tokenize('_')
    println z
    it[1] = z[0]+':' + z[1]

    if(!aMap[it[0]].contains(it[1])) {
      aMap[it[0]] << it[1]
    }

    cList << it[1]
    println "${++t}/${aFileContent.size()}"
    }
  }
}

if(additionalAssocation && additionalAssocation != 'false') {
  new File(additionalAssocation).splitEachLine('\t') {
    def z = it[1].tokenize('/').last().tokenize('_')
    it[1] = z[0] + ':' + z[1]

    if(!aMap[it[0]]) { aMap[it[0]] = [] }

    if(!aMap[it[0]].contains(it[1])) {
      aMap[it[0]] << it[1]
    }

    cList << it[1]
  }
}

println 'writing the annotation file now'
def sWriter = new BufferedWriter(new FileWriter('annot.tsv'))
def oo = ''
def y = 0
aMap.each { a, b ->
  println "(${++y}/${aMap.size()})"
  if(!b.any{ it.indexOf('txt') != -1}) {
    sWriter.write('http://reality.rehab/ptvis/' + a + '\t' + b.join(';') + '\n')
  }
}
sWriter.flush()
sWriter.close()
println 'done'

cList = cList.unique()

println cList

URIFactory factory = URIFactoryMemory.getSingleton()

def ontoFile = oFile
URI graph_uri = factory.getURI("http://gloo/")
['LOO', 'DOID', 'BFO', 'RO', 'ECO', 'UPHENO', 'GENO', 'UBERON', 'NCBITaxon', 'SYMP', 'HP', 'CL', 'TRANS', 'CHEBI', 'SO'].each {
  factory.loadNamespacePrefix(it, factory.getURI("http://"+it + '/').toString())
}

slib.graph.model.graph.G graph = new GraphMemory(graph_uri)

//gConf.addGDataConf(dataConf)
//gConf.addGAction(actionRerootConf)

def annot = 'annot.tsv'
GDataConf aConf = new GDataConf(GFormat.TSV_ANNOT, annot)
GraphLoaderGeneric.populate(aConf, graph)

GDataConf dataConf = new GDataConf(GFormat.RDF_XML, oFile)
GraphLoaderGeneric.populate(dataConf, graph)

// Add virtual root for 3 subontologies__________________________________
URI virtualRoot = factory.getURI("http://gloo/virtualRoot")
graph.addV(virtualRoot)

GAction rooting = new GAction(GActionType.REROOTING)
rooting.addParameter("root_uri", virtualRoot.stringValue())
GraphActionExecutor.applyAction(factory, rooting, graph)

/*
def graphURI = factory.getURI('http://DOID/')
factory.loadNamespacePrefix("DOID", graphURI.toString());

G graph = new GraphMemory(graphURI)

def dataConf = new GDataConf(GFormat.RDF_XML, ontoFile)
def actionRerootConf = new GAction(GActionType.REROOTING)
//actionRerootConf.addParameter("root_uri", "HP:0000118"); // phenotypic abnormality
actionRerootConf.addParameter("root_uri", "DOID:4"); // phenotypic abnormality

def gConf = new GraphConf()
//gConf.addGDataConf(dataConf)
//gConf.addGAction(actionRerootConf)
def annot = 'annot.tsv'
gConf.addGDataConf(new GDataConf(GFormat.TSV_ANNOT, annot));

GraphLoaderGeneric.load(gConf, graph)
*/

println graph.toString()

def roots = new ValidatorDAG().getTaxonomicRoots(graph)
println roots

def icConf = new IC_Conf_Corpus(SMConstants.FLAG_IC_ANNOT_RESNIK_1995)
//def icConf = new IC_Conf_Topo(SMConstants.FLAG_ICI_ZHOU_2008)
//def icConf = new IC_Conf_Topo(SMConstants.FLAG_ICI_SANCHEZ_2011)
def smConfPairwise = new SMconf(SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_RESNIK_1995, icConf)
//def smConfPairwise = new SMconf(SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_LIN_1998, icConf)
//def smConfGroupwise = new SMconf(SMConstants.FLAG_SIM_GROUPWISE_AVERAGE, icConf)
def smConfGroupwise = new SMconf(SMConstants.FLAG_SIM_GROUPWISE_BMA, icConf)
// FLAG_SIM_GROUPWISE_AVERAGE_NORMALIZED_GOSIM

//def smConfPairwise = new SMconf(SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_JIANG_CONRATH_1997_NORM , icConf)


def out = []
def z = 0

def outWriter = new BufferedWriter(new FileWriter(aFile + '_sim_matrix_fast5_2500_'+oFile.tokenize('/').last()+'.lst'), 1024 * 1024 * 1024)
if(rmneg) {
outWriter = new BufferedWriter(new FileWriter(aFile + '_sim_matrix_fast5_noneg.lst'), 1024 * 1024 * 1024)
}

def engine = new SM_Engine(graph)

cList = cList.unique()


def getURIfromTerm = { term ->
    term = term.tokenize(':')

    //println term
    //return factory.getURI("http://${term[0]}/" + term[1])
    //println ''
    def start = 'purl.obolibrary.org/obo/'
    if(term[0] == 'LOO') {
      start = 'reality.rehab/ontologies/'
    }
    return factory.getURI("http://$start" + term[0] + '_' + term[1])
}

def bestMatches = [:]
aMap.each { g1, u1 ->
  def bestGroup = []
  def bestCount = 0

  aMap.each { g2, u2 ->
    def matchingDiagnoses = patient_visit_diagnoses[g1].findAll { v ->
      patient_visit_diagnoses[g2].contains(v)
    }.size()

    if(matchingDiagnoses == bestCount) {
      bestGroup << g2
    }
    if(matchingDiagnoses < bestCount) {
      bestGroup = [ g2 ]
      bestCount = matchingDiagnoses 
    }
  }

  bestMatches[g1] = bestGroup
}


def rrs=[]
def aps=[]
aMap.each { g1, u1 ->
  println "(${++z}/${aMap.size()})"

  def aList = []
  aMap.each { g2, u2 ->
    def match = patient_visit_diagnoses[g1][0] == patient_visit_diagnoses[g2][0]
    def match2 = patient_visit_diagnoses[g1].any { v -> patient_visit_diagnoses[g2].contains(v) }
    def match3 = patient_visit_diagnoses[g2].any { patient_visit_diagnoses[g1][0] == it }
    def match4 = bestMatches[g1].contains(g2)

try {
    aList << [g1,g2,engine.compare(smConfGroupwise, smConfPairwise,
                                    u1.collect { 
                                      getURIfromTerm(it)
                                     }.findAll { graph.containsVertex(it) }.toSet(), 
                                    u2.collect { 
                                      getURIfromTerm(it)
                                    }.findAll { graph.containsVertex(it) }.toSet())
                                    
                                   // u1.collect { factory.getURI(it) }.toSet(), 
                                    //u2.collect { factory.getURI(it) }.toSet())
                                    ,match,match2,match3,match4]
                                  } catch(e) { aList << [g1,g2,0,match,match2,match3,match4] }
  }
  aList = aList.toSorted { it[2] }.reverse() //[0..10]

  def matchRank = 0
  aList.eachWithIndex { it, i ->
    if(matchRank > 0) { return; }
    if(patient_visit_diagnoses[it[0]][0] == patient_visit_diagnoses[it[1]][0]) {
      matchRank = i+1 
    }
  }
  rrs << 1/matchRank //

  aList.eachWithIndex { it, i -> 
    outWriter.write(it.join(',') + ',' + (i+1) + '\n')
  }

  aList = aList[0..10]

  def ps = []  // precisions
  // for each patient it[1] similar to it[0], see if it[1] contains the diagnosis from it[0] in position 0
  aList.eachWithIndex { it, i ->
    if(patient_visit_diagnoses[it[1]][0] == patient_visit_diagnoses[it[0]][0]) {
      ps << (ps.size()+1) / (i+1) //
    }
  }
  def ap
  try {
    ap = ps.sum() / ps.size() //
  } catch(e) { ap = 0 }
  println ap
  aps << ap  
}

def mrr = rrs.sum() / rrs.size() //
def map = aps.sum() / aps.size() //
println 'mrr: ' + mrr
println 'map: ' + map

outWriter.flush()
outWriter.close()

/*
self-implementation of BMA (not used in experiment, but it can be faster for larger amounts of patients)
def cSim = [:] // class similarity
cList.eachWithIndex { u1, i ->
  println "${i}/${cList.size()}"
  cSim[u1] = [:]
  cList.each { u2 ->
    cSim[u1][u2] = engine.compare(smConfPairwise, 
      factory.getURI('http://purl.obolibrary.org/obo/'+u1.replace(':','_')),
      factory.getURI('http://purl.obolibrary.org/obo/'+u2.replace(':','_')))
  }
}

println 'doing bma'

GParsPool.withPool(45) { p ->  // im a monster
aMap.eachParallel { g1, u1 ->
  println "(${++z}/${aMap.size()})"
  def aList = []
  aMap.each { g2, u2 ->
    aList << [
      g2,
      ((u1.inject(0) { sum, uri ->
        sum += u2.collect { uri2 -> cSim[uri][uri2] }.max()
      } + u2.inject(0) { sum, uri ->
        sum += u1.collect { uri2 -> cSim[uri][uri2] }.max()
      }) / (u1.size() + u2.size())) // /
    ]
  }

  aList.toSorted { it[1] }.reverse().eachWithIndex { i, it ->
    outWriter.write(g1 + ',' + it[0] + ',' + (i+1) + ','+ it[1]+ ',' + match + ','+match2+'\n')
  }
}
}*/

outWriter.flush()
outWriter.close()

//new File(aFile + '_sim_matrix.lst').text = out.join('\n')

//sum += u2.inject(0) { max, uri2 -> (cSim[uri][uri2] > max) ? cSim[uri][uri2] : max }
//sum += u1.inject(0) { max, uri2 -> (cSim[uri][uri2] > max) ? cSim[uri][uri2] : max }
