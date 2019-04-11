import DtAlignment.DriftTube as DriftTube
import DtAlignment.DtModule as DtModule
import DtAlignment.utils
import ROOT

def display_event(track,dt_modules):
    """ Draw a simple event in x-z-projection with positions of all the tubes in the
    detector, the drift circles as well as the reconstructed track.
    Note that this does not take RPCs into account.
    
    Parameters
    ----------
    track:
        genfit Track object containing the fit as well as the points used for fitting
    
    dt_modules:
        dictionary containing all the DtModule objects of which the detector is built
    """
    canvas = ROOT.TCanvas("Simple DT display")
    
    #Draw all the tubes
    tubes = []
    for key in dt_modules.keys():
        for tube in dt_modules[key].get_tubes():
            tubes.append(ROOT.TEllipse(tube._position[0],tube._position[2],1.85))
    
    for tube in tubes:
        tube.Draw()
        
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
        hits.append(ROOT.TEllipse(tube._position[0],tube._position[2],drift_radius))
        
    for hit in hits:
        hit.Draw()
        
    #TODO draw track line
    canvas.SaveAs("Test_disp.pdf")
    
    