# -*- coding: utf-8 -*-
"""Doc
GFF3, GTF
col 9s
col 1-8
col 9


- class Father # class 装在文件
- class Sun # 解析文件细节
"""
import gzip

class Annotation:
    """Father classs of GTF or GFF3.
    Loader of files.
    """
    REQUIRED_COLUMNS = [
        "seq_name",
        "source",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]

    def __init__(self,file):
        # loading file
        self.__file = file
        self._lines = None
        self.__load_file__()
        pass

    def __load_file__(self):
        self._lines = open(self.__file,"rt") if ".gz" not in self.__file else gzip.open(self.__file,"rt")


    pass


class Gff3(Annotation):
    """Gff3 class.
    Parse Gff3 files.
    """
    def to_gtf(self,file):
        # convert gff3 to gtf
        print("Convert Gff3 obj to Gtf obj...")
        f = open(file, "wt") if ".gz" not in file else gzip.open(file,"wt")
        # print(self._lines.readline())

        counter_line = 0

        for line in self._lines:
            line = line.strip()
            counter_line+=1

            if line:
                if line.startswith("#")
                    print(f'Skip #annotation line {counter_line}')
                    continue
                else: #core
                    try:
                        seq_name,source,feature,start,end,score,strand,frame,attribute = line.split("\t")
                    except Exception as e:
                        raise gtftools.ParsingError()

                    attribute ={
                        k.rstrip(): v.strip()
                        for k,v in
                        [i.split("=") for i in attribute.split(";")]
                    }
                    _dt ={
                        k:v
                        for k,v in zip(
                            self.REQUIRED_COLUMNS,
                            (seq_name,source,feature,start,end,\
                             score,strand,frame,attribute)
                        )
                    }
                    pass


        f.close()

class Gtf(Annotation):
    """Gff3 class.
    Parse Gff3 files.
    """
    pass
