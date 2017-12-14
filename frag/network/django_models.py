from django.db import models

class Edge(models.Model):
    """
    Django model to hold Edge information
    """
    # Now parse the relevant CSV file and bulk add the data
    exclude_smiles = models.CharField(max_length=200)
    exclude_type = models.CharField(max_length=200)
    rebuilt_smiles = models.CharField(max_length=200)
    rebuilt_ring_smiles = models.CharField(max_length=200)
    rebuilt_type = models.CharField(max_length=200)
    excluded_ring_smiles = models.CharField(max_length=200)
    # The nodes
    node_from = models.ForeignKey(Nodes)
    node_to = models.ForeignKey(Nodes)
    # Define this column as unique
    unique_together = (("exclude_smiles", "rebuilt_smiles"),)



class Nodes(models.Model):
    """
    Django model to hold Node information
    """
    smiles = models.CharField(unique=True,max_length=200)
    heavy_atom_count = models.IntegerField()
    ring_atom_count = models.IntegerField()
    ring_smiles = models.CharField(max_length=200)
    # Then add the annotations
    price = models.IntegerField()
    mol_type = models.CharField(max_length=200)

