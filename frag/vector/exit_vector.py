from rdkit import Chem
import numpy as np
from numpy import linalg


def get_best_fit_plane(pts, weights=None):
    """

    :param pts:
    :param weights:
    :return:
    """
    if weights is None:
        wSum = len(pts)
        origin = np.sum(pts, 0)
    origin /= wSum
    sums = np.zeros((3, 3), np.double)
    for pt in pts:
        dp = pt - origin
        for i in range(3):
            sums[i, i] += dp[i] * dp[i]
            for j in range(i + 1, 3):
                sums[i, j] += dp[i] * dp[j]
                sums[j, i] += dp[i] * dp[j]
    sums /= wSum
    vals, vects = linalg.eigh(sums)
    order = np.argsort(vals)
    normal = vects[:, order[0]]
    plane = np.zeros((4, ), np.double)
    plane[:3] = normal
    plane[3] = -1 * normal.dot(origin)
    return plane

def plane_best_fit_calc(mol, confId=-1):
    """

    :param mol:
    :param confId:
    :return:
    """
    conf = mol.GetConformer(confId)
    if not conf.Is3D():
        return 0
    pts = np.array([list(conf.GetAtomPosition(x)) for x in range(mol.GetNumAtoms())])
    plane = get_best_fit_plane(pts)
    denom = np.dot(plane[:3], plane[:3])
    denom = denom**0.5
    # add up the distance from the plane for each point:
    res = 0.0
    for pt in pts:
        res += np.abs(pt.dot(plane[:3]) + plane[3])
        res /= denom
        res /= len(pts)
    return res

def plane_best_fit_exit_vector(scaffold):
    """

    :param scaffold:
    :return:
    """

    # Get PBF plane for murcko scaffold only
    confId = -1
    conf = scaffold.GetConformer(confId)
    if not conf.Is3D():
        print('This mol is not 3D - all PBFev angles will be 0 degrees')
        return [0]
    pts = np.array([list(conf.GetAtomPosition(i))  # Get atom coordinates if not Xenon
                    for i in xrange(scaffold.GetNumAtoms()) if scaffold.GetAtomByIdx(i).element != "Xe"])
    # Plane is xyz vector with a c intercept adjustment
    plane = get_best_fit_plane(pts)

    # Where [#0] matches exit vector SMILES [*]
    patt = Chem.MolFromSmarts('[#0]-[Xe]')
    matches = scaffold.GetSubstructMatches(patt)
    if len(matches) == 0:
        return None

    # Calculate angles between exit vectors and the murcko plane of best fit
    exitVectors = np.zeros(len(matches))
    denom = np.dot(plane[:3], plane[:3])
    denom = denom**0.5
    for n, match in enumerate(matches):
        evCoords = conf.GetAtomPosition(match[0])
        anchorCoords = conf.GetAtomPosition(match[1])
        v = np.array(((evCoords[0]-anchorCoords[0]),
                      (evCoords[1]-anchorCoords[1]),
                      (evCoords[2]-anchorCoords[2])))
        angle = np.arcsin((np.dot(v, plane[:3])) /
                          ((denom)*((np.dot(v, v))**0.5)))
        angle = np.abs(np.degrees(angle))
        exitVectors[n] = angle
    return exitVectors

