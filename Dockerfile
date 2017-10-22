FROM informaticsmatters/rdkit:Release_2017_03_3
RUN pip install tqdm nose neo4j-driver
ADD . /usr/local/fragalysis
RUN pip install /usr/local/fragalysis
