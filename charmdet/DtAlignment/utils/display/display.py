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
    
    #TODO not yet drawn
    _draw_detector(canvas, dt_modules)
        
    #Draw hits
    hits = []
    for point in track.getPointsWithMeasurement():
        raw = point.getRawMeasurement()
        det_id = raw.getDetId()
        parse_result = DtAlignment.utils.parse_det_id(det_id)
        drift_radius = raw.getRawHitCoords()[6]
        module = dt_modules[parse_result['module']]
        #TODO horribly slow
        for j in range(len(module.get_tubes())):
            tube = module.get_tubes()[j]
            if tube._ID == det_id:
                break
        tube = module.get_tubes()[j]
        hit_circle = ROOT.TEllipse(tube._position[2],tube._position[0],drift_radius)
        hit_circle.SetFillColor(ROOT.kRed)
        hits.append(hit_circle)
        
    for hit in hits:
        hit.Draw()
        
    trajectory = _create_trackline(track)
    trajectory.Draw()
    canvas.SaveAs("Test_disp.pdf")
    
def _draw_detector(canvas,dt_modules):
    """ Draws the whole detector into the passed canvas object.
    
    Parameters
    ----------
    canvas: ROOT.TCanvas
        TCanvas or TPad in that the detector should be drawn
        
    dt_modules: dict
        Dictionary with module names as keys and DtModule objects contained. Everything in this dict will be drawn
    """
    #canvas.cd()
    tubes = []
    for key in dt_modules.keys():
        for tube in dt_modules[key].get_tubes():
            tubes.append(ROOT.TEllipse(tube._position[2],tube._position[0],1.815))
    
    for tube in tubes:
        tube.Draw()
    
def _draw_trackline(canvas,track):
    """ Draws a line for the track trajectory onto the passed TCanvas or TPad object.
    This draws __segments__ of the track between each pair of two consecutive hits.
    
    Parameters
    ----------
    canvas: ROOT.TCanvas
        TCanvas or TPad in that the detector should be drawn
    
    track: genfit.Track
        Track object that should be drawn
    """
    #canvas.cd()
    n_points = track.getNumPointsWithMeasurement()
    trajectory_points = [0] * n_points
    trajectory_momenta = [0] * n_points
    for i in range(n_points):
        fitted_state = track.getFittedState(i)
        trajectory_momenta[i] = fitted_state.getMom()
        trajectory_points[i] = fitted_state.getPos()
    
       
    x_coords = np.empty(n_points,dtype=np.float64)
    z_coords = np.empty(n_points,dtype=np.float64)
    for i in range(n_points):
        x_coords[i] = trajectory_points[i][0]
        z_coords[i] = trajectory_points[i][2]
    polyline = ROOT.TPolyLine(n_points,z_coords,x_coords)
    polyline.Draw()
    
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
    
    return ROOT.TPolyLine(n_points,z_coords,x_coords)
        
    