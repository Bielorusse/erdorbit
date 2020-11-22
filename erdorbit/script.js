/*
script.js for Erdorbit
author: thibaut voirand
*/

/*Description of the script*************************************************************************
The purpose of this script is to compute the positions of a
spacecraft on a given orbit in cartesian coordinates, and to draw
them.
***************************************************************************************************/

/*Declaration of variables**************************************************************************
***************************************************************************************************/

// declaring  calculation variables
var pos = [[]]; // positions vectors array
var pos2D = [[]]; // projected positions vectors array
var MU_PLANET = 398600.4418; // planet gravitational parameter (km3/s2)
var ROT_PLANET = 360 / 86164.1004; // planet rotation velocity (deg/s)
var dur; // simulation duration
var dt; // simulation step-time
var numberOfSteps; // simulation number of steps
var a; // orbit semi-major axis
var e; // orbit eccenctricity
var i; // orbit inclination
var RAAN = 0; // orbit right ascension of the ascending node
var om; // orbit argument of periapsis

// declaring drawing related variables
var DRAWING_SIZE_FACTOR = 9 / 10;
var canvas = document.getElementById('Canvas');
var canvasContext = canvas.getContext('2d');
canvasContext.lineWidth = 0.3;
var timer; // timer variable for animation
var alpha = 0 * Math.PI / 180; // plane angle around intersection of itself and OXY plane
var beta = 0 * Math.PI / 180; // plane angle around OZ
var delta = 0 * Math.PI / 180; // rotation of projected figure around normal to projection plane
var xTranslation = 0; // translation along x-axis of object from center of drawing
var yTranslation = 0; // translation along y-axis of object from center of drawing

/*Math- geometry- and canvas-related functions******************************************************
***************************************************************************************************/

function fromOrbitalToCartesianCoordinates(a, e, i, RAAN, om, t, MU_PLANET) {
  /*
	Converting from orbital parameters to cartesian coordinates
	- Inputs:
			a       	semi-major axis (km)
			e       	eccentricity (-)
			i       	inclination (deg)
			RAAN    	right ascension of the ascending node (deg)
			om      	argument of periapsis (deg)
			t       	time spent since passage at periapsis (s)
			MU_PLANET	gravitational parameter of the planet	(km3/s2)
	- Outputs:
			pos   	  position vector of the orbiting object (km)
	*/

  // converting angles to radians
  i = i * Math.PI / 180;
  RAAN = RAAN * Math.PI / 180;
  om = om * Math.PI / 180;

  // computing mean anomaly
  var n = Math.sqrt(MU_PLANET / Math.pow(a, 3));
  var M = n * t;

  // computing eccentric anomaly
  var E = [M];
  var j;
  for (j = 0; j < 100; j++) {
    E[j + 1] =
      E[j] + (M - E[j] + e * Math.sin(E[j])) / (1 - e * Math.cos(E[j]));
    if (Math.abs(E[j + 1] - E[j]) < 1e-8) {
      var E = E[j + 1];
      break;
    }
  }

  // computing true anomaly
  var nu =
    2 *
    Math.atan2(
      Math.sqrt(1 + e) * Math.sin(E / 2),
      Math.sqrt(1 - e) * Math.cos(E / 2)
    ) %
    (Math.PI * 2);

  // computing radius
  var r = a * (1 - Math.pow(e, 2)) / (1 + e * Math.cos(nu));

  // computing pos vector
  var pos = [
    r *
      (Math.cos(om + nu) * Math.cos(RAAN) -
        Math.sin(om + nu) * Math.sin(RAAN) * Math.cos(i)),
    r *
      (Math.cos(om + nu) * Math.sin(RAAN) +
        Math.sin(om + nu) * Math.cos(RAAN) * Math.cos(i)),
    r * (Math.sin(om + nu) * Math.sin(i))
  ];

  // returning result
  return pos;
}

function fromJ2000ToECEF(posJ2000, t, ROT_PLANET) {
  /*
	Converting coordinates from standard J2000 to Earth-centered Earth-fixed reference frame
  - Inputs:
      posJ2000    position vector in rotating frame (km)
      t 		 	    time (s)
      ROT_PLANET	planet rotational velocity (deg/s)
  - Outpus:
      posECEF	    position vector in planet fixed ref. frame (km)
	*/

  ROT_PLANET = ROT_PLANET * Math.PI / 180;

  var posECEF = [
    Math.cos(ROT_PLANET * t) * posJ2000[0] +
      Math.sin(ROT_PLANET * t) * posJ2000[1],
    -Math.sin(ROT_PLANET * t) * posJ2000[0] +
      Math.cos(ROT_PLANET * t) * posJ2000[1],
    posJ2000[2]
  ];

  return posECEF;
}

function getExtremum(array) {
  /*
	Getting extremum of a 2D array
	- Inputs:
	   array     2D array
	- Outputs:
	   extremum  extremum value of 2D array
	*/

  var j; // creating count variable
  var k; // creating count variable

  var extremum = array[0][0];

  for (j = 0; j < array.length; j++) {
    for (k = 0; k < array[j].length; k++) {
      if (extremum < Math.abs(array[j][k])) {
        extremum = Math.abs(array[j][k]);
      }
    }
  }

  return extremum;
}

function multiplyArrayByScalar(inputArray, scal) {
  /*
	Multiplying each component of a 2D array by a scalar
	- Inputs:
      inputArray  input array
      scal        input scalar
	- Outputs:
		  outputArray output array
	*/

  var outputArray = new Array(inputArray.length);

  var j; // creating count variable
  var k; // creating count variable
  for (j = 0; j < inputArray.length; j++) {
    outputArray[j] = new Array(inputArray[j].length);
    for (k = 0; k < inputArray[j].length; k++) {
      outputArray[j][k] = inputArray[j][k] * scal;
    }
  }
  return outputArray;
}

function planarProjection(inputArray, alpha, beta) {
  /*
	Projecting an array of 3D coordinates on a plane
  - Inputs:
	   inputArray 	 array of 3D coordinates
		 beta	         rotation of the projection plane around Oz
		 alpha         rotation of the p.p. around intersection of itself and Oxy plane
		 delta         rotation of the projected vector around the normal to the projection plane
	- Outputs:
	   outputArray   2D projection of the input array on the plane
	*/

  var outputArray = new Array(inputArray.length);

  var j;
  for (j = 0; j < inputArray.length; j++) {
    outputArray[j] = new Array(2);
    outputArray[j][0] =
      Math.cos(beta) * inputArray[j][0] - Math.sin(beta) * inputArray[j][1];
    outputArray[j][1] =
      -Math.sin(beta) * Math.sin(alpha) * inputArray[j][0] -
      Math.cos(beta) * Math.sin(alpha) * inputArray[j][1] +
      Math.cos(alpha) * inputArray[j][2];
  }

  return outputArray;
}

function computePositions(
  a,
  e,
  i,
  RAAN,
  om,
  MU_PLANET,
  ROT_PLANET,
  numberOfSteps,
  dt
) {
  /*
  This function computes the positions of the orbiting object in a ECEF reference frame
  - Inputs:
    - orbital parameters:
      a           semi-major axis
      e           eccentricity
      i           inclination
      RAAN        right ascension of the ascending node
      om          argument of periapsis
    - planetary constants:
      MU_PLANET   planet gravitational parameter
      ROT_PLANET  planet rotational velocity
    numberOfSteps number of steps
    dt            step time
  Ouputs:
    positionsArray array of positions
  */

  var positionsArray = new Array(numberOfSteps);

  // computing positions for each simulation step
  var j; // creating count variable
  for (j = 0; j < numberOfSteps; j++) {
    positionsArray[j] = fromOrbitalToCartesianCoordinates(
      a,
      e,
      i,
      RAAN,
      om,
      j * dt,
      MU_PLANET
    );
    positionsArray[j] = fromJ2000ToECEF(positionsArray[j], j * dt, ROT_PLANET);
  }

  return positionsArray;
}

function adaptCoordinatesToCanvasFrame(inputCoord, canvasWidth, canvasHeight) {
  /*
  The reference frame of the HTML5 canvas is centered in the upper left corner, with the X axis
  pointing to the right, and the Y axis pointing down.
  This function transforms 2D coordinates so that the drawing is centered at the center of the
  canvas, with the Y axis pointing up.
  - Input:
    inputCoord   array of 2D coordinates expressed in a "classical" cartesian reference frame
    canvasWidth  canvas width
    canvasHeight canvas height
  - Output:
    outputCoord  array of 2D coordinates, adapted to the HTML5 canvas' reference frame
  */

  var outputCoord = new Array(inputCoord.length);

  var j; // creating count variable
  for (j = 0; j < inputCoord.length; j++) {
    outputCoord[j] = new Array(2);
    outputCoord[j][0] = inputCoord[j][0] + canvasWidth / 2;
    outputCoord[j][1] = -inputCoord[j][1] + canvasHeight / 2;
  }

  return outputCoord;
}

function translateCoordinates(inputCoord, dx, dy) {
  /*
  This function translates the positions in a 2D reference frame
  - Inputs:
      inputCoord  input coordinates
      dx          translation along x axis
      dy          translation along y axis
  - Output:
      outputCoord output coordinates
  */

  var outputCoord = new Array(inputCoord.length);

  var j;
  for (j = 0; j < inputCoord.length; j++) {
    outputCoord[j] = new Array(2);
    outputCoord[j][0] = inputCoord[j][0] + dx;
    outputCoord[j][1] = inputCoord[j][1] + dy;
  }

  return outputCoord;
}

function rotateCoordinates(inputCoord, delta) {
  /*
  This functions rotates 2D coordinates around the center of the reference frame from a delta angle
  - Input:
      inputCoord  input 2D coordinates
      delta       rotation angle
  - Output:
      outputCoord output 2D coordinates
  */

  var outputCoord = new Array(inputCoord.length);

  var j;
  for (j = 0; j < inputCoord.length; j++) {
    outputCoord[j] = new Array(2);
    outputCoord[j][0] =
      Math.cos(delta) * inputCoord[j][0] + Math.sin(delta) * inputCoord[j][1];
    outputCoord[j][1] =
      Math.cos(delta) * inputCoord[j][1] - Math.sin(delta) * inputCoord[j][0];
  }

  return outputCoord;
}

function resizeDrawingToFitCanvas(inputCoord, canvasHeight) {
  /*
  This function resizes the drawing to fit the canvas size
  - Input:
      inputCoord    input coordinates
      canvasHeight  canvas height
  - Output:
      outputCoord   output coordinates
  */

  var drawingSize = getExtremum(inputCoord) * 2;

  var outputCoord = multiplyArrayByScalar(
    inputCoord,
    canvasHeight / drawingSize * DRAWING_SIZE_FACTOR
  );

  return outputCoord;
}

function drawOrbit(
  pos,
  alpha,
  beta,
  delta,
  xTranslation,
  yTranslation,
  canvasContext,
  canvasWidth,
  canvasHeight
) {
  /*
  This function draws the trajectory of an orbiting object on the canvas
  This function calls several previously defined functions
  - Input:
      pos   array of 3D cartesian positions the orbiting object
      alpha         projection plane angle
      beta          projection plane angle
      delta         drawing rotation angle
      xTranslation  drawing translation along x-axis
      yTranslation  drawing translation along y-axis
      canvasContext canvas context
      canvasWidth   canvas width
      canvasHeight  canvas height
  */

  canvasContext.clearRect(0, 0, canvasWidth, canvasHeight);

  pos2D = planarProjection(pos, alpha, beta);

  pos2D = rotateCoordinates(pos2D, delta);

  pos2D = translateCoordinates(pos2D, xTranslation, yTranslation);

  pos2D = adaptCoordinatesToCanvasFrame(pos2D, canvasWidth, canvasHeight);

  // Drawing line
  canvasContext.beginPath();
  var j; // creating count variable
  canvasContext.moveTo(pos2D[0][0], pos2D[0][1]);
  for (j = 1; j < pos2D.length; j++) {
    canvasContext.lineTo(pos2D[j][0], pos2D[j][1]);
  }
  canvasContext.stroke();
}

function getInputParameters() {
  /*
  This function assigns parameters values entered by the user in the input forms to the
  calculation variables
  */

  // simulation parameters

  // simulation duration (s)
  dur = parseFloat(document.getElementById('DurationInput').value * 86400);
  // step time (s) and number of steps (-)
  if (document.getElementById('StepTimeRadioButton').checked) {
    dt = parseFloat(document.getElementById('StepTimeInput').value);
    numberOfSteps = dur / dt;
  } else {
    numberOfSteps = parseInt(
      document.getElementById('NumberOfStepsInput').value,
      10
    );
    dt = dur / numberOfSteps;
  }

  // orbital parameters

  // semi-major axis (km)
  a = parseFloat(document.getElementById('SemiMajorAxisInput').value);
  // eccentricity (-)
  e = parseFloat(document.getElementById('EccentricityInput').value);
  // inclination (deg)
  i = parseFloat(document.getElementById('InclinationInput').value);
  // argument of periapsis (deg)
  om = parseFloat(document.getElementById('ArgumentPeriapsisInput').value);
}

/*User-Interaction-related functions****************************************************************
***************************************************************************************************/

document.getElementById('DrawButton').onclick = function drawButton() {
  /*
  This function gets the user-entered input parameters, computes the corresponding orbiting object
  positions, resizes the positions coordinates to fit the canvas size, and draws the trajectory on
  the canvas
	*/

  getInputParameters();

  pos = computePositions(
    a,
    e,
    i,
    RAAN,
    om,
    MU_PLANET,
    ROT_PLANET,
    numberOfSteps,
    dt
  );

  pos = resizeDrawingToFitCanvas(pos, canvas.height);

  drawOrbit(
    pos,
    alpha,
    beta,
    delta,
    xTranslation,
    yTranslation,
    canvasContext,
    canvas.width,
    canvas.height
  );
};

document.getElementById('ClearButton').onclick = function clearButton() {
  /*
  This functions clears the canvas and sets the plane angles and drawing translations to zero
  */

  canvasContext.clearRect(0, 0, canvas.width, canvas.height);
  alpha = 0;
  beta = 0;
  delta = 0;
  xTranslation = 0;
  yTranslation = 0;
};

document.getElementById('CenterButton').onclick = function centerButton() {
  /*
  This function recenters the drawing by setting the translation and angle variables to zero
  It then draws the orbit again
  */

  if (document.getElementById('TranslateCheckbox').checked) {
    // resetting object shift from center to zero
    xTranslation = 0;
    yTranslation = 0;
  } else {
    // resetting projection plane direction angles to zero
    alpha = 0;
    beta = 0;
    delta = 0;
  }

  drawOrbit(
    pos,
    alpha,
    beta,
    delta,
    xTranslation,
    yTranslation,
    canvasContext,
    canvas.width,
    canvas.height
  );
};

document.getElementById('UpButton').onclick = function upButton() {
  /*
  This function either:
    - executes a positive rotation of the projection plane around the intersection of
      itself and of the OXY plane, by increasing
    - or translates the drawing along the y-axis
  */

  if (document.getElementById('TranslateCheckbox').checked) {
    yTranslation = yTranslation - 5; // shifting object from center
  } else {
    alpha = alpha + 5 * Math.PI / 180; // rotating projection plane
  }

  drawOrbit(
    pos,
    alpha,
    beta,
    delta,
    xTranslation,
    yTranslation,
    canvasContext,
    canvas.width,
    canvas.height
  );
};

document.getElementById('DownButton').onclick = function downButton() {
  /*
  This function either:
    - executes a negative rotation of the projection plane around the intersection of
      itself and of the OXY plane, by increasing
    - or translates the drawing along the y-axis
  */

  if (document.getElementById('TranslateCheckbox').checked) {
    yTranslation = yTranslation + 5; // shifting object from center
  } else {
    alpha = alpha - 5 * Math.PI / 180; // rotating projection plane
  }

  drawOrbit(
    pos,
    alpha,
    beta,
    delta,
    xTranslation,
    yTranslation,
    canvasContext,
    canvas.width,
    canvas.height
  );
};

document.getElementById('LeftButton').onclick = function leftButton() {
  /*
  This function either:
    - executes a positive rotation of the projection plane around axis OZ, by increasing
    - or translates the drawing along the x-axis
  */

  if (document.getElementById('TranslateCheckbox').checked) {
    xTranslation = xTranslation + 5; // shifting object from center
  } else {
    beta = beta + 5 * Math.PI / 180; // rotating projection plane
  }

  drawOrbit(
    pos,
    alpha,
    beta,
    delta,
    xTranslation,
    yTranslation,
    canvasContext,
    canvas.width,
    canvas.height
  );
};

document.getElementById('RightButton').onclick = function rightButton() {
  /*
  This function either:
    - executes a negative rotation of the projection plane around axis OZ, by increasing
    - or translates the drawing along the x-axis
  */

  if (document.getElementById('TranslateCheckbox').checked) {
    xTranslation = xTranslation - 5; // shifting object from center
  } else {
    beta = beta - 5 * Math.PI / 180; // rotating projection plane
  }

  drawOrbit(
    pos,
    alpha,
    beta,
    delta,
    xTranslation,
    yTranslation,
    canvasContext,
    canvas.width,
    canvas.height
  );
};

document.getElementById(
  'ClockwiseButton'
).onclick = function clockwiseButton() {
  /*
  This function executes a positive rotation of the coordinates around the normal to the projection
  plane
  */

  delta = delta + 5 * Math.PI / 180; // rotating projection plane

  drawOrbit(
    pos,
    alpha,
    beta,
    delta,
    xTranslation,
    yTranslation,
    canvasContext,
    canvas.width,
    canvas.height
  );
};

document.getElementById(
  'CounterclockwiseButton'
).onclick = function counterclockwiseButton() {
  /*
  This function executes a negative rotation of the coordinates around the normal to the projection
  plane
  */

  delta = delta - 5 * Math.PI / 180; // rotating projection plane

  drawOrbit(
    pos,
    alpha,
    beta,
    delta,
    xTranslation,
    yTranslation,
    canvasContext,
    canvas.width,
    canvas.height
  );
};

document.getElementById('ZoomInButton').onclick = function zoomInButton() {
  /*
  This function executes a zoom in
  */

  // increasing values of positions array
  pos = multiplyArrayByScalar(pos, 1.1);

  drawOrbit(
    pos,
    alpha,
    beta,
    delta,
    xTranslation,
    yTranslation,
    canvasContext,
    canvas.width,
    canvas.height
  );
};

document.getElementById('ZoomOutButton').onclick = function zoomOutButton() {
  /*
  This function executes a zoom out
  */

  // decreasing values of positions array
  pos = multiplyArrayByScalar(pos, 0.9);

  drawOrbit(
    pos,
    alpha,
    beta,
    delta,
    xTranslation,
    yTranslation,
    canvasContext,
    canvas.width,
    canvas.height
  );
};

document.getElementById(
  'AdaptSizeButton'
).onclick = function adaptSizeButton() {
  /*
  This function resizes the drawing to fit the canvas size
  */

  pos = resizeDrawingToFitCanvas(pos, canvas.height);

  drawOrbit(
    pos,
    alpha,
    beta,
    delta,
    xTranslation,
    yTranslation,
    canvasContext,
    canvas.width,
    canvas.height
  );
};

document.getElementById(
  'PlayAnimationButton'
).onclick = function playAnimationButton() {
  /*
  This function displays the orbiting object's trajectory step by step in an animation
  */

  // clearing the canvas
  canvasContext.clearRect(0, 0, canvas.width, canvas.height);

  // drawing a new point every millisecond
  function animate() {
    var j = 0;
    timer = setInterval(function() {
      canvasContext.beginPath();
      canvasContext.moveTo(pos2D[j][0], pos2D[j][1]);
      canvasContext.lineTo(pos2D[j + 1][0], pos2D[j + 1][1]);
      canvasContext.stroke();
      j++;
    }, 1);
  }
  animate();
};

document.getElementById(
  'StopAnimationButton'
).onclick = function stopAnimationButton() {
  /*
  This function stops the animation
  */

  clearTimeout(timer);
};

document.getElementById('DownloadButton').onclick = function downloadButton() {
  /*
	In this function, a new canvas with a different resolution is created
	The trajectory is plotted again on this new canvas, and it is downloaded
	*/

  if (document.getElementById('DownloadResolution').value == 600) {
    // creating a 600x600px canvas
    var canvasToDownload = document.createElement('canvas');
    canvasToDownload.width = 600;
    canvasToDownload.height = 600;
  } else if (document.getElementById('DownloadResolution').value == 1000) {
    // creating a 1000x1000px canvas
    var canvasToDownload = document.createElement('canvas');
    canvasToDownload.width = 1000;
    canvasToDownload.height = 1000;
  } else if (document.getElementById('DownloadResolution').value == 2000) {
    // creating a 2000x2000px canvas
    var canvasToDownload = document.createElement('canvas');
    canvasToDownload.width = 2000;
    canvasToDownload.height = 2000;
  } else if (document.getElementById('DownloadResolution').value == 10000) {
    // creating a 10000x10000px canvas
    var canvasToDownload = document.createElement('canvas');
    canvasToDownload.width = 10000;
    canvasToDownload.height = 10000;
  }

  //creating canvas context
  canvasToDownloadContext = canvasToDownload.getContext('2d');
  canvasToDownloadContext.lineWidth =
    0.3 * canvasToDownload.height / canvas.height; // setting line width

  // adjusting image size to new canvas resolution
  posCanvasToDownload = multiplyArrayByScalar(
    pos,
    canvasToDownload.height / canvas.height
  );

  drawOrbit(
    posCanvasToDownload,
    alpha,
    beta,
    delta,
    xTranslation * canvasToDownload.width / canvas.width,
    yTranslation * canvasToDownload.height / canvas.height,
    canvasToDownloadContext,
    canvasToDownload.width,
    canvasToDownload.height
  );

  // creating a (fictional) link to download the file
  var downloadLink = document.createElement('a');
  downloadLink.download = 'Erdorbit.png'; // file name

  //storing canvas data in a blob object for downloading
  canvasToDownload.toBlob(function(blob) {
    downloadLink.href = URL.createObjectURL(blob); //creating data URL
    downloadLink.click(); //
  }, 'image/png');
};

document.getElementById('UploadButton').onclick = function uploadButton() {
  /*
	In this function, the parameters entered by the user are stored
	in a variable, and sent to the server along with the canvas data
	through a php script.
	The php script stores the parameters in a .csv file, and stores
	the canvas data in a .png file in a gallery folder
	*/

  // converting drawing projection plane angles to degrees
  var alpha_deg = alpha * 180 / Math.PI;
  var beta_deg = beta * 180 / Math.PI;
  var delta_deg = delta * 180 / Math.PI;

  // creating string containing drawing parameters
  var drawing_param =
    String(a) +
    ',' +
    String(e) +
    ',' +
    String(i) +
    ',' +
    String(om) +
    ',' +
    String(dur / 86400) +
    ',' + // simulation duration is stored in days
    String(numberOfSteps);

  // getting username
  var username = document.getElementById('UsernameInput').value;
  var useremail = document.getElementById('UserEmailInput').value;

  // getting URL of canvas data
  var canvasDataURL = canvas.toDataURL();

  // ajax request to send data to server
  $.ajax({
    type: 'POST', // type of ajax request best fit for large amount of data
    url: 'upload.php', // URL to send request to
    data: {
      imgBase64: canvasDataURL, // sending canvas data encrypted in base 64
      username: username,
      useremail: useremail,
      drawing_param: drawing_param
    }
  }).done(function() {
    alert('upload complete'); // displaying success to the user
    console.log('upload complete'); // displaying success in the console
  });
};

/*Gallery-related ajax request**********************************************************************
***************************************************************************************************/

// displaying all png images from the gallery directory to the html page
$.ajax({
  // retrieving database info from server with ajax request
  type: 'GET',
  url: 'uploads_database.csv', // location of database
  success: function(data) {
    // transferring drawings info from server to javascript array
    var csv = []; // array that will receive drawings info
    var data_lines = data.split(/\r\n|\n/); // split lines from csv database
    var headers = data_lines[0].split(','); // extract headers (first line)

    for (var i = 1; i < data_lines.length; i++) {
      var data_cell = data_lines[i].split(','); // split cells from csv
      if (data_cell.length == headers.length) {
        var single_line = [];
        for (var j = 0; j < headers.length; j++) {
          single_line.push(data_cell[j]);
        }
        csv.push(single_line);
      }
    }

    // append drawings info to html file
    for (j = 0; j < csv.length; j++) {
      if (csv[j][2]) {
        // in case a username was given
        $('#GalleryDiv').append(
          "<img src='gallery/" +
          'Bild' +
          csv[j][0] +
          '_' +
          csv[j][1] +
          '_' +
          csv[j][4] +
          ".png'><br>" + // append image link
          "<div class='Tooltip'>" +
          'Bild' +
          csv[j][0] +
          '_' +
          csv[j][1] + // image label text
          "<span class='TooltipText'>" + // tooltip image info
            '<u>Uploaded by:</u><br>' +
            csv[j][2] +
            '<br>' +
            '<u>Parameters:</u><br>' +
            'a(km):' +
            csv[j][5] +
            '<br>' +
            'e(-):' +
            csv[j][6] +
            '<br>' +
            'i(deg):' +
            csv[j][7] +
            '<br>' +
            'om(deg):' +
            csv[j][8] +
            '<br>' +
            'dur(days):' +
            csv[j][9] +
            '<br>' +
            'step(-):' +
            csv[j][10] +
            '<br>' +
            '</span></div><br><br>'
        );
      } else {
        // in cas no username was given
        $('#GalleryDiv').append(
          "<img src='gallery/" +
          'Bild' +
          csv[j][0] +
          '_' +
          csv[j][1] +
          '_' +
          csv[j][4] +
          ".png'><br>" + // append image link
          "<div class='Tooltip'>" +
          'Bild' +
          csv[j][0] +
          '_' +
          csv[j][1] + // image label text
          "<span class='TooltipText'>" + // tooltip image info
            '<u>Parameters:</u><br>' +
            'a(km):' +
            csv[j][5] +
            '<br>' +
            'e(-):' +
            csv[j][6] +
            '<br>' +
            'i(deg):' +
            csv[j][7] +
            '<br>' +
            'om(deg):' +
            csv[j][8] +
            '<br>' +
            'dur(days):' +
            csv[j][9] +
            '<br>' +
            'step(-):' +
            csv[j][10] +
            '<br>' +
            '</span></div><br><br>'
        );
      }
    }
  }
});
