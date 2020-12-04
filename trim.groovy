def mapped_doids = new File('./komenti_input.txt').text.split('\n')

def trimmed_unexp = []
new File('./disease_vocab.txt').splitEachLine('\t') {
  if(mapped_doids.contains(it[1])) { trimmed_unexp << it.join('\t') }
}
new File('trim_disease_vocab.txt').text = trimmed_unexp.join('\n')

