#!/usr/bin/env groovy
@Grab(group='commons-cli', module='commons-cli', version='1.4')
@Grab(group='org.apache.commons', module='commons-lang3', version='3.4')
@Grab(group='edu.stanford.nlp', module='stanford-corenlp', version='3.7.0')
@Grab(group='edu.stanford.nlp', module='stanford-corenlp', version='3.7.0', classifier='models')
@Grab(group='edu.stanford.nlp', module='stanford-parser', version='3.7.0')

import java.util.concurrent.*
import java.util.concurrent.atomic.*
import groovyx.gpars.*
import org.codehaus.gpars.*

import edu.stanford.nlp.pipeline.*
import edu.stanford.nlp.ling.*
import edu.stanford.nlp.semgraph.*

def props = new Properties()
props.put("annotators", "tokenize, ssplit, pos, lemma")
StanfordCoreNLP pipeline = new StanfordCoreNLP(props)

def files = []
new File('./texts').eachFile { files << it }

def i = 0
GParsPool.withPool(85) { p ->  // he he he
  files.eachParallel{ e ->
    def text = e.text
    text = text.replaceAll('\n\n', '.\n')
    text = text.replaceAll('\u2022', '.  ')
    text = text.replaceAll('â€“', '. ')
    text = text.replaceAll('-', '.  ')
    text = text.replaceAll('– ', '. ')
    text = text.replaceAll('\\s+', ' ')
    text = text.replaceAll(', \\?', '. ?')
    text = text.replaceAll('\\.', '. ')

    def aDocument = new Annotation(text.toLowerCase())
    pipeline.annotate(aDocument)

    def newText = ''
    aDocument.get(CoreAnnotations.SentencesAnnotation.class).each { sentence ->
      newText += sentence.get(CoreAnnotations.TokensAnnotation.class).collect {
        it.lemma()
      }.join(' ') + ' '
    }

    newText = newText.replaceAll(' ,', ',')
    newText = newText.replaceAll(' \\.', '.')

/*    if(newText.indexOf('discharge diagnosis') != -1) {
      newText = newText.replaceAll('discharge diagnosis*', '')
*/
      new File('new_texts/' + e.getName()).text = newText
 //   }

    println "${++i}/${files.size()}"
  }
}

