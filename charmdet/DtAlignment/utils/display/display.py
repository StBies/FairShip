import DtAlignment.DriftTube as DriftTube
import DtAlignment.DtModule as DtModule
import DtAlignment.utils
import ROOT
import numpy as np

def display_event(track,dt_modules):
    """ Draw a simple event in x-z-projection with positions of all the tubes in the
    detector, the drift circles as well as the reconstructed track.
    Note that this does not take RPCs into account.
    
    Parameters
    ----------
    track: genfit.Track
        genfit Track object containing the fit as well as the points used for fitting
    
    dt_modules: dict
        dictionary containing all the DtModule objects of which the detector is built
    """
    canvas = ROOT.TCanvas("Simple DT display")
    canvas.Range(-20,-120,800,120)
    
    tubes = []
    for key in dt_modules.keys():
        for tube in dt_modules[key].get_tubes():
            tubes.append(ROOT.TEllipse(tube._position[2],tube._position[0],1.815))
    
    for tube in tubes:
        tube.Draw("SAME")
        
    hits = _create_driftcircles(track, dt_modules) 
    
    for hit in hits:
        hit.Draw()
        
    trajectory = _create_trackline(track)
    trajectory.Draw()
    canvas.SaveAs("Test_disp.pdf")
  
#TODO doesn't work probably due to wrong pointer to canvas  
def _draw_detector(canvas,dt_modules):
    """ Draws the whole detector into the passed canvas object.
    
    Parameters
    ----------
    canvas: ROOT.TCanvas
        TCanvas or TPad in that the detector should be drawn
        
    dt_modules: dict
        Dictionary with module names as keys and DtModule objects contained. Everything in this dict will be drawn
    """
    canvas.cd()
    tubes = []
    for key in dt_modules.keys():
        for tube in dt_modules[key].get_tubes():
            tubes.append(ROOT.TEllipse(tube._position[2],tube._position[0],1.815))
    
    for tube in tubes:
        tube.Draw("SAME")
    
    
def _create_trackline(track):
    """ Create a TPolyLine from the hits of a given track object.
    
    Parameters
    ----------
    track: genfit.Track
        Track object that should be drawn
        
    Returns
    -------
    ROOT.TPolyLine
        Object containing the line information that is drawable
    """
    n_points = track.getNumPointsWithMeasurement()
    
    x_coords = np.empty(n_points,dtype=np.float64)
    z_coords = np.empty(n_points,dtype=np.float64)
    
    for i in range(n_points):
        fitted_state = track.getFittedState(i)
        x_coords[i] = fitted_state.getPos()[0]
        z_coords[i] = fitted_state.getPos()[2]
    
    polyline = ROOT.TPolyLine(n_points,z_coords,x_coords)
    polyline.SetLineColor(ROOT.kMagenta)
    
    return polyline

def _create_driftcircles(track,dt_modules):
    """ Creates and returns a list of drift circles that can then be drawn to any canvas.
    
    Parameters
    ----------
    track: genfit.Track
        Track for that the hits should be drawn
        
    dt_modules: dict
        dictionary containing the DtModule objects. The coordinates of the hit wires are read from here.
        
    Returns
    -------
    List of ROOT.TEllipse objects
    """
    hits = []
    for point in track.getPointsWithMeasurement():
        raw = point.getRawMeasurement()
        det_id = raw.getDetId()
        parse_result = DtAlignment.utils.parse_det_id(det_id)
        drift_radius = raw.getRawHitCoords()[6]
        module = dt_modules[parse_result['module']]
        #TODO horribly slow, could read det coords from hit.
        for j in range(len(module.get_tubes())):
            tube = module.get_tubes()[j]
            if tube._ID == det_id:
                break
        tube = module.get_tubes()[j]
        hit_circle = ROOT.TEllipse(tube._position[2],tube._position[0],drift_radius)
        hit_circle.SetFillColor(ROOT.kRed)
        hit_circle.SetLineWidth(1)
        hits.append(hit_circle)
        
    return hits
        
    