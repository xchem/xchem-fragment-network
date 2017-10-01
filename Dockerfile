FROM informaticsmatters/rdkit
RUN pip install tqdm
ADD . /usr/local/fragalysis
RUN pip install /usr/local/fragalysis
