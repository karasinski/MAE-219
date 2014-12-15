#### import the simple module from the paraview
from paraview.simple import *
import sys
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
cavityClippedfoam = GetActiveSource()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1054, 790]

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')
pLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 5e-17, 0.865003, 0.865003, 0.865003, 1e-16, 0.705882, 0.0156863, 0.14902]
pLUT.ScalarRangeInitialized = 1.0

# show data in view
cavityClippedfoamDisplay = Show(cavityClippedfoam, renderView1)
# trace defaults for the display properties.
cavityClippedfoamDisplay.ColorArrayName = ['POINTS', 'p']
cavityClippedfoamDisplay.LookupTable = pLUT
cavityClippedfoamDisplay.ScalarOpacityUnitDistance = 1.0844426982393176
cavityClippedfoamDisplay.SelectInputVectors = ['POINTS', 'U']
cavityClippedfoamDisplay.WriteLog = ''

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
cavityClippedfoamDisplay.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')
pPWF.Points = [0.0, 0.0, 0.5, 0.0, 1e-16, 1.0, 0.5, 0.0]
pPWF.ScalarRangeInitialized = 1

# change representation type
cavityClippedfoamDisplay.SetRepresentationType('Surface LIC')

# set scalar coloring
ColorBy(cavityClippedfoamDisplay, ('POINTS', 'streamFunction'))

# rescale color and/or opacity maps used to include current data range
cavityClippedfoamDisplay.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
cavityClippedfoamDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'streamFunction'
streamFunctionLUT = GetColorTransferFunction('streamFunction')
streamFunctionLUT.RGBPoints = [-2.037569999694824, 0.231373, 0.298039, 0.752941, -0.9894677493721247, 0.865003, 0.865003, 0.865003, 0.058634500950574875, 0.705882, 0.0156863, 0.14902]
streamFunctionLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'streamFunction'
streamFunctionPWF = GetOpacityTransferFunction('streamFunction')
streamFunctionPWF.Points = [-2.037569999694824, 0.0, 0.5, 0.0, 0.058634500950574875, 1.0, 0.5, 0.0]
streamFunctionPWF.ScalarRangeInitialized = 1

# Properties modified on renderView1
renderView1.Background = [1.0, 1.0, 1.0]

# Properties modified on renderView1
renderView1.OrientationAxesVisibility = 0

# Properties modified on cavityClippedfoamDisplay
cavityClippedfoamDisplay.ColorMode = 'Multiply'

# hide color bar/color legend
cavityClippedfoamDisplay.SetScalarBarVisibility(renderView1, False)

# Properties modified on cavityClippedfoamDisplay
cavityClippedfoamDisplay.EnhanceContrast = 'LIC and Color'

#change interaction mode for render view
renderView1.InteractionMode = '2D'

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [-6.546626850823488e-06, 0.024711792637900302, 0.29411509449255]
renderView1.CameraFocalPoint = [-6.546626850823488e-06, 0.024711792637900302, 0.004999999888241291]
renderView1.CameraParallelScale = 0.0768830580193085

#### uncomment the following to render all views
SaveScreenshot('/Users/localmin/code/MAE-219/cavityClipped/test_image.png', magnification=1, quality=100, view=renderView1)
sys.exit()
