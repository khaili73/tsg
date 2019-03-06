#!/usr/bin/python

data = open("./data.fastq", "r")
filtered = open("./filtered.fastq", "w")
#head = [next(data) for x in xrange(12)]

N_seq_ind = []
i = 0
record = []
to_remove = False
total_score = 0

for line in data:
   i += 1
   record.append(line)
   if i % 4 == 2:
      if "N" in line:
         to_remove = True
   if i % 4 == 0:
      for q in line:
         total_score += ord(q) - 33
      avg_score = total_score / len(line)
      total_score = 0
      if not to_remove and avg_score >= 20:
         for r in record:
            filtered.write(r)
      to_remove = False
      record = []

filtered.close()
data.close()
