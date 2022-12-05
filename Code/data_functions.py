#########################
# File with different function definitions, used in converting HyDAMO data to a D-HyDAMO model

def get_crosssection_culvert_AGV(shape: int = 1, height: float = None, 
                                width: float = None, closed:int = 1):
    '''
    Function defines a cross section based on different parameters that interact with the shape of the culvert.

    '''
    shapedict = {1:"circle",3: "rectangle"} 
    shape_str = shapedict[shape]

    # Include the diameter when the culvert is a circle
    if shape == 1: diameter = width
    else: diameter = None

    crosssection = {"shape"    : shape_str,
                    "diameter" : diameter,
                    "height"   : height,
                    "width"    : width,
                    "closed"   : closed,
                    }
    return crosssection

