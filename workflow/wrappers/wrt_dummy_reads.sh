zcat $1 | awk --posix '{ if (NR % 4 == 0) { sub(/.*/,"CCCCCCCC")} else if (NR % 2 == 0) { sub(/.*/,"NNNNNNNN")}; print}' | gzip > $2
