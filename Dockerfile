FROM informaticsmatters/rdkit
RUN pip install tqdm nose neo4j-driver
ADD . /usr/local/fragalysis
RUN pip install /usr/local/fragalysis
