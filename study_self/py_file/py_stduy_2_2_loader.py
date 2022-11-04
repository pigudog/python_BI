"""
FASTQ 读取
FASTA 读取
"""

import gzip


def load_fastx(file) -> list:
  """
  :param str file: path of input
  :return: a list for all reads
    i.e. [['header','seq','info','quality'],[]]
  :rtype: list
  """
  f = open(file,"rt") if not ".gz" in file else gzip.open(file,"rt")
  # FASTQ @ ; FASTA >. 判断文件类型
  symbol = f.read(1)
  print(symbol)
  f.close()
  # 全部载入内存就好，当文件很小的时候
  f = open(file, "rt") if not ".gz" in file else gzip.open(file, "rt")
  if symbol == "@":
    # FASTQ
    print("is fastq!")
    raw_info = [i.rstrip() for i in f.readlines()]
    # print(raw_info)
    ls = []
    read = []
    for line in raw_info:
      if line.startswith("@"):
        # header line!
        n = len(read)
        if n == 0: ## read == []
          # first circle
          read.append(line)
        elif n == 4: ##read=['header','seq','info','quality']
          ls.append(read)
          read = []
          read.append(line)
        else:
          f.close()
          raise ValueError("The file may be incomplete")
      else:
        # not header line!
        read.append(line)

  elif symbol == ">":
    # FASTA
    # print('is fasta!')
    ls = []
    read = []
    seq = ''
    raw_info = [i.rstrip() for i in f.readlines()]
    # print(raw_info)

    for line in raw_info:
      if line.startswith('>'):
        # 读取 header！
        n = len(read)

        if n == 0:
          # 第一次循环
          read.append(line)
        elif n == 1:
          # 已经有一个 header 了！现在缺 seq
          read.append(seq)  # add seq line
          ls.append(read)
          read = []  # 重置 read 这个 list
          read.append(line)
          seq = ''  # 重置 seq 这个 str
        else:
          f.close()
          raise ValueError('The file may be incomplete!')
      else:
        # 读取并添加 seq
        seq += line

    read.append(seq)
    ls.append(read)
  else:
    raise ValueError("Input line one must starts with @ for FQ or > for FA")
  return ls
  f.close()

  # 当文件很大的时候，用生成器函数产生一个迭代器进行文件的读取
def load_fastx_generato(file):
  """
  :param str file: path of input
  :return: a list for all reads
    i.e. print(next(obj)) -> ['header', 'seq', 'info', 'quality']
  :rtype: list
  """
  f = open(file, "rt") if not ".gz" in file else gzip.open(file, "rt")
  # FASTQ @ ; FASTA >. 判断文件类型
  symbol = f.read(1)
  print(symbol)
  f.close()
  # 全部载入内存就好，当文件很小的时候
  f = open(file, "rt") if not ".gz" in file else gzip.open(file, "rt")
  if symbol == "@":
    # FASTQ
    print("is fastq!")
    # print(raw_info)
    ls = []
    read = []
    line = f.readline().rstrip()
    while True:
      if not line: # readline 后期会全是空格的返回
        break
      else:
        if line.startswith('@'):
          n = len(read)
          if n == 0:
            # read == []
            # 第一次循环开始
            read.append(line)  # add header! 注意使用后我们需要调用下一行
            line = f.readline().rstrip()
          elif n == 4:
            # read == ['header', 'seq', 'info', 'quality']
            # ls.append(read)
            yield read
            read = []
            read.append(line)
            line = f.readline().rstrip()
          elif n == 3:
            # TODO
            """
            @@IIEIBCE>IC<IBIIIIEAIEIEB<IDECCD6 # line! 期望它是 header！现在它是 quality
            [
                '@Beta12AdemL1C001R00100001768/1', 
                'ATCCCCGTATCTTCACCCCACCACAAACTATTAG', 
                '+',
            ]'@@IIEIBCE>IC<IBIIIIEAIEIEB<IDECCD6'

            """
            read.append(line)
            line = f.readline().rstrip()

            if not line.startswith('@'):
              raise ValueError('The file may be incomplete!')
          else:
            f.close()
            raise ValueError('The file may be incomplete!')

          # header line!
          # read.append(line)  # add header !
        else:
          # not header line!
          read.append(line)
          line = f.readline().rstrip()

  elif symbol == ">":
    # FASTA
    # print('is fasta!')
    ls = []
    read = []
    seq = ''
    raw_info = [i.rstrip() for i in f.readlines()]
    # print(raw_info)

    for line in raw_info:
      if line.startswith('>'):
        # 读取 header！
        n = len(read)

        if n == 0:
          # 第一次循环
          read.append(line)
        elif n == 1:
          # 已经有一个 header 了！现在缺 seq
          read.append(seq)  # add seq line
          ls.append(read)
          read = []  # 重置 read 这个 list
          read.append(line)
          seq = ''  # 重置 seq 这个 str
        else:
          f.close()
          raise ValueError('The file may be incomplete!')
      else:
        # 读取并添加 seq
        seq += line

    read.append(seq)
    ls.append(read)


if __name__ == '__main__':
  FQ_TEST = "fake_fq.fastq"
  FA_TEST = "fake_fa.fasta"
  fq = load_fastx(file=FQ_TEST)
  fa = load_fastx(file=FA_TEST)
  print(fq,fa)
  pass
