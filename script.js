
/*
script.js for Erdorbit
author: thibaut voirand
*/



/*Description of the script*****************************************

The purpose of this script is to compute the positions of a
spacecraft on a given orbit in cartesian coordinates, and to draw
them.

The first part of the script contains:
	- mathematical functions for the calculations
	- declaration of variables
	- function which plots the positions on a HTML5 canvas

The second part of the script contains functions that allow the user
to manipulate the drawing, such as:
	- zoom in and out
	- rotate the projection plane (and thus, the image of the trajectory)
	- clear, recenter the canvas
	- play an animation of the trajectory

The third part of the script contains functions that allow the user
to store the drawing (download to computer or upload to gallery), as well
as a request to display all images from the gallery on the page

*******************************************************************/





/*First part of the script******************************************

The first part of the script contains mathematical functions for the
calculations, declaration of variables, and the "draw" function
which plots the positions on a HTML5 canvas.

Mathematical functions:
	- orb2car		converts orbital parameters to cartesian coordinates
	- rot_frame		converts coordinates to rotating, Earth-fixed frame
	- getmax		returns greatest value from a 2D array
	- list_prod		multiply each value from an array by a scalar
	- proj			projects spatial coordinates on a plane: 3D to 2D

Variables:
	Calculation:
	- pos			array containing positions' 3D coordinates
	- pos_2d		array containing positions' 2D coordinates
	- mu_planet		planet's gravitational parameter
	- rot_planet	planet's rotation velocity
	- dur			duration of simulation
	- dt			simulation step-time
	- step_nb		number of steps
	- a				orbit's semi-major axis
	- e				orbit's eccentricity
	- i				orbit's inclination
	- RAAN			orbit's right ascension of the ascending node
	- om			orbit's argument of periapsis
	Drawing:
	- size			drawing size, scope of the canvas
	- canvas		HTML5 canvas
	- ctx			canvas' context
	- cheight		canvas' height
	- timer			timer for animation
	- alpha			rotation angle of projection plane along
					intersection of itself and OXY plane
	- beta			rotation angle of projection plane along axis OZ
	- delta			rotation of projected coordinates around
					normal to projection plane
	- dx			translation of the coordinates along X
	- dy			translation of the coordinates along Y

*******************************************************************/



// Additional functions


// converting from orbital parameters to cartesian coordinate
function orb2car(a, e, i, RAAN, om, t, mu_planet){

	/*

	Converting from orbital parameters to cartesian coordinates

		Inputs:
			a       	semi-major axis (km)
			e       	eccentricity (-)
			i       	inclination (deg)
			RAAN    	right ascension of the ascending node (deg)
			om      	argument of periapsis (deg)
			t       	time spent since passage at periapsis (s)
			mu_planet	gravitational parameter of the planet
						(km3/s2)

		Outputs:
			pos   	position vector of the satellite (km)

	*/

    // converting angles to radians
    i = i * Math.PI/180 ;
    RAAN = RAAN * Math.PI/180 ;
    om = om * Math.PI/180 ;

	// computing mean anomaly
	var n = Math.sqrt(mu_planet/ Math.pow(a,3)) ;
	var M = n * t ;

	// computing eccentric anomaly
	var E = [M] ;
	var j ;
	for (j=0; j<100; j++ ) {
		E[j+1] = E[j] + (M - E[j] + e * Math.sin(E[j])) /
			(1 - e * Math.cos(E[j])) ;
		if (Math.abs(E[j+1] - E[j]) < 1e-8){
			var E = E[j+1] ;
			break
		}
	}

	// computing true anomaly
	var nu = (2 * Math.atan2( Math.sqrt( (1+e)/(1-e) ) *
		Math.tan(E/2) , 1 ) ) % 360 ;

	// computing radius
	var r = a * (1 - e * Math.cos(E) );

	// computing focal parameter
	var p = a * (1 - Math.pow(e,2));

	// computing specific angular momentum
	var h = Math.sqrt( mu_planet * p ) ;

	// computing pos vector
	var pos = [r * (Math.cos(om+nu) * Math.cos(RAAN) -
					Math.sin(om+nu) * Math.sin(RAAN) * Math.cos(i)),
				r * (Math.cos(om+nu) * Math.sin(RAAN) +
				Math.sin(om+nu) * Math.cos(RAAN) * Math.cos(i)),
				r * (Math.sin(om+nu) * Math.sin(i))] ;

	// returning result
    return pos;

}

// converting from rotating reference frame to planet fixed  r. f.
function rot_frame(pos,t,rot_planet){

	/*

	Converting from rotating to planet fixed reference frame

	Inputs:
		pos			position vector in rotating frame (km)
		t 		 	time (s)
		rot_planet	planet rotational velocity (deg/s)

	Outpus:
		pos_rot	=	position vector in planet fixed ref. frame (km)

	*/

	rot_planet = rot_planet * Math.PI/180 ;

	var pos_rot = [Math.cos(rot_planet*t) * pos[0] +
				Math.sin(rot_planet*t) * pos[1] ,
				- Math.sin(rot_planet*t) * pos[0] +
				Math.cos(rot_planet*t) * pos[1],
				pos[2]] ;

	return pos_rot ;

}

// getting max value of a 2D array
function getmax(array){
	/*
	Getting max value of a 2D array

	Inputs:
		array 	= 2D array

	Outputs:
		max 	= max value of 2D array

	*/

	var j ; // creating count variable
	var k ; // creating count variable

	var max = array[0][0] ;

	for (j=0 ; j<array.length ; j++){
		for (k=0 ; k<array[j].length ; k++){
			if (max < array[j][k]){
				max = array[j][k] ;
			}
		}
	}

	return max;
}

// mutiplying each cell of an array by a scalar
function list_prod(list,scal){

	/*

	Multiplying each component of a list of lists by a scalar

	Inputs:
		list 	= input list of lists
		scal	= input scalar

	Outputs:
		list 	= output vector

	*/

	var j ; // creating count variable
	var k ; // creating count variable
	for (j=0 ; j<list.length ; j++){
		for (k=0 ; k<list[j].length ; k++){
			list[j][k] = list[j][k] * scal ;
		}
	}
	return list ;
}

// projecting a 3D vector on a plane
function proj(vec, alpha, beta, delta){

	/*
	Projecting a 3D vector on a plane

	Inputs:
		vec 	= 3D vector
		beta	= rotation of the projection plane around Oz
		alpha    = rotation of the p.p. around intersection
					of itself and Oxy plane
		delta   = rotation of the projected vector around the normal
				  to the projection plane

	Outputs:
		vec 	= 2D projection of the input vector on the plane

	*/


	vec = [Math.cos(beta) * vec[0] -
					Math.sin(beta) * vec[1],
					- Math.sin(beta) * Math.sin(alpha) * vec[0]
					- Math.cos(beta) * Math.sin(alpha) * vec[1]
					+ Math.cos(alpha) * vec[2]];

	vec = [Math.cos(delta) * vec[0] + Math.sin(delta) * vec[1],
			Math.cos(delta) * vec[1] - Math.sin(delta) * vec[0]];



	return vec ;

}


// declaring  calculation variables
var pos = [[]] ; // positions vectors array
var pos_2d = [[]] ; // projected positions vectors array
var mu_planet ; // planet gravitational parameter
var rot_planet ; // planet rotation velocity
var dur ; // simulation duration
var dt ; // simulation step-time
var step_nb ; // simulation number of step
var a ; // orbit semi-major axis
var e ; // orbit eccenctricity
var i ; // orbit inclination
var RAAN ; // orbit right ascension of the ascending node
var om ; // orbit argument of periapsis

// Setting up the drawing
var size ; // drawing size
var canvas = document.getElementById("canvas");
var ctx = canvas.getContext("2d");
ctx.lineWidth = 0.3 ;
var cheight = canvas.height ;
var timer ; // declaring timer variable for animation

// declaring projection plane direction angles
var alpha = 0 * Math.PI/180 ;// plane angle around intersection of
// itself and OXY plane
var beta = 0 * Math.PI/180 ;// plane angle around OZ
var delta = 0 * Math.PI/180 ; // rotation of the projected figure
// around the normal to the projection plane

// declaring translation of object from center of drawing
var dx = 0 ;
var dy = 0 ;


// drawing orbit when clicking the "draw" button
document.getElementById("bt_draw").onclick = function (){

	/*
	Drawing orbit

	In this function:
		- the parameters (orbit, simulation, and planet) are entered
		- the positions are computed and saved in the pos[[]] array
		- the positions are adapted (size, 2D projection, reference frame)
		and drawn on the canvas


	*/


	// Setting up the problem


	// planet parameters

	// pl. gravitational parameter (km3/s2)
	mu_planet = 399000 ;
	// pl. rotational velocity (deg/s)
	rot_planet = 0.00417 ;

	// simulation parameters

	// simulation duration (s)
	dur = document.getElementById("input_dur").value*86400 ;
	// step time (s) and number of steps (-)
	if (document.getElementById("radio_step_time").checked){
		dt = document.getElementById("step_time").value ;
		step_nb = dur/dt ;
	} else {
		step_nb = document.getElementById("step_nb").value ;
		dt = dur/step_nb ;
	}

	// orbital parameters

	// semi-major axis (km)
	a = document.getElementById("input_a").value ;
	// eccentricity (-)
	e = document.getElementById("input_e").value ;
	// inclination (deg)
	i = document.getElementById("input_i").value ;
	// right ascension of the ascending node (deg)
	RAAN = 0 ;
	// argument of periapsis (deg)
	om = document.getElementById("input_om").value ;


	// getting positions for each simulation step

	var j ; // creating count variable
	for (j=0 ; j<step_nb ; j++){
		pos[j] = orb2car(a,e,i,RAAN,om,j*dt,mu_planet);
		pos[j] = rot_frame(pos[j],j*dt,rot_planet);
	}

	// getting size of the drawing
	size = getmax(pos)*2;

	// resizing positions to fit the graph
	pos = list_prod(pos, cheight/size*9/10) ;

	// Projecting points on 2D plane
	var j ; // creating count variable
	for(j=0 ; j<pos.length ; j++){
		pos_2d[j] = (proj(pos[j], alpha, beta, delta)) ;
	}

	// converting positions reference frame to fit the HTML canvas
	var j ; // creating count variable
	for (j=0 ; j<pos_2d.length ; j++){
		pos_2d[j][0] = pos_2d[j][0]+cheight/2 +dx;
		pos_2d[j][1] = -pos_2d[j][1]+cheight/2 -dy;
	}

	// clearing the canvas
	ctx.clearRect(0, 0, cheight, cheight) ;

	// Drawing line
	ctx.beginPath() ;
	var j ; // creating count variable
	ctx.moveTo(pos_2d[0][0], pos_2d[0][1]);
	for (j=1 ; j<pos_2d.length ; j++){
		ctx.lineTo(pos_2d[j][0], pos_2d[j][1]);
	}
	ctx.stroke();

}





/*Second part of the script*****************************************

This part of the script contains functions that allow the user to
manipulate the drawing.

Function:
	- clearing canvas
	- recentering drawing
	- positive rotation of projection plane around intersection of
	itself and OXY ("up")
	- negative rotation of projection plane around intersection of
	itself and OXY ("down")
	- positive rotation of projection plane around axis OZ ("left")
	- negative rotation of projection plane around axis OZ ("right")
	- positive rotation of coordinates along normal to projection plane
	("clockwise")
	- negative rotation of coordinates along normal to projection plane
	("counterclockwise")
	- zoom in
	- zoom out
	- adapt size of the drawing to size of the canvas
	- play an animation of the trajectory
	- stop animation of the trajectory


*******************************************************************/


// clearing the canvas when clicking the "clear" button
document.getElementById("bt_clear").onclick = function(){

	ctx.clearRect(0, 0, cheight, cheight) ;
	alpha = 0 ;
	beta = 0 ;
	delta = 0 ;
	dx = 0 ;
	dy = 0 ;

}

// recentering the drawing when clicking the "center" button
document.getElementById("bt_center").onclick = function() {

	if(document.getElementById("checkbox_translate").checked){
		// resetting object shift from center to zero
		dx = 0;
		dy = 0;
	} else {
		// resetting projection plane direction angles to zero
		alpha = 0 ;
		beta = 0 ;
		delta = 0 ;
	}

	// clearing the canvas
	ctx.clearRect(0, 0, cheight, cheight) ;

	// Projecting points on 2D plane
	var j ; // creating count variable
	for(j=0 ; j<pos.length ; j++){
		pos_2d[j] = (proj(pos[j], alpha, beta, delta)) ;
	}

	// converting positions to fit the HTML canvas frame
	var j ; // creating count variable
	for (j=0 ; j<pos_2d.length ; j++){
		pos_2d[j][0] = pos_2d[j][0]+cheight/2 +dx;
		pos_2d[j][1] = -pos_2d[j][1]+cheight/2 -dy;
	}

	// Drawing line
	ctx.beginPath() ;
	var j ; // creating count variable
	ctx.moveTo(pos_2d[0][0], pos_2d[0][1]);
	for (j=1 ; j<pos_2d.length ; j++){
		ctx.lineTo(pos_2d[j][0], pos_2d[j][1]);
	}
	ctx.stroke();

}

// "up": positive rotation of projection plane around intersection of itself and OXY
document.getElementById("bt_up").onclick = function() {

	if(document.getElementById("checkbox_translate").checked){
		dy = dy - 5 ; // shifting object from center
	} else {
		alpha = alpha + 5 * Math.PI/180 // rotating projection plane
	}

	// clearing the canvas
	ctx.clearRect(0, 0, cheight, cheight) ;

	// Projecting points on 2D plane
	var j ; // creating count variable
	for(j=0 ; j<pos.length ; j++){
		pos_2d[j] = (proj(pos[j], alpha, beta, delta)) ;
	}

	// converting positions to fit the HTML canvas frame
	var j ; // creating count variable
	for (j=0 ; j<pos_2d.length ; j++){
		pos_2d[j][0] = pos_2d[j][0]+cheight/2 +dx;
		pos_2d[j][1] = -pos_2d[j][1]+cheight/2 -dy;
	}

	// Drawing line
	ctx.beginPath() ;
	var j ; // creating count variable
	ctx.moveTo(pos_2d[0][0], pos_2d[0][1]);
	for (j=1 ; j<pos_2d.length ; j++){
		ctx.lineTo(pos_2d[j][0], pos_2d[j][1]);
	}
	ctx.stroke();


}

// "down": negative rotation of projection plane around intersection of itself and OXY
document.getElementById("bt_down").onclick = function() {

	if(document.getElementById("checkbox_translate").checked){
		dy = dy + 5 ; // shifting object from center
	} else {
		alpha = alpha - 5 * Math.PI/180 // rotating projection plane
	}

	// clearing the canvas
	ctx.clearRect(0, 0, cheight, cheight) ;

	// Projecting points on 2D plane
	var j ; // creating count variable
	for(j=0 ; j<pos.length ; j++){
		pos_2d[j] = (proj(pos[j], alpha, beta, delta)) ;
	}

	// converting positions to fit the HTML canvas frame
	var j ; // creating count variable
	for (j=0 ; j<pos_2d.length ; j++){
		pos_2d[j][0] = pos_2d[j][0]+cheight/2 +dx;
		pos_2d[j][1] = -pos_2d[j][1]+cheight/2 -dy;
	}

	// Drawing line
	ctx.beginPath() ;
	var j ; // creating count variable
	ctx.moveTo(pos_2d[0][0], pos_2d[0][1]);
	for (j=1 ; j<pos_2d.length ; j++){
		ctx.lineTo(pos_2d[j][0], pos_2d[j][1]);
	}
	ctx.stroke();

}

// "left": positive rotation of projection plane around axis OZ
document.getElementById("bt_left").onclick = function() {

	if(document.getElementById("checkbox_translate").checked){
		dx = dx + 5 ; // shifting object from center
	} else {
		beta = beta + 5 * Math.PI/180 // rotating projection plane
	}

	// clearing the canvas
	ctx.clearRect(0, 0, cheight, cheight) ;

	// Projecting points on 2D plane
	var j ; // creating count variable
	for(j=0 ; j<pos.length ; j++){
		pos_2d[j] = (proj(pos[j], alpha, beta, delta)) ;
	}

	// converting positions to fit the HTML canvas frame
	var j ; // creating count variable
	for (j=0 ; j<pos_2d.length ; j++){
		pos_2d[j][0] = pos_2d[j][0]+cheight/2 +dx;
		pos_2d[j][1] = -pos_2d[j][1]+cheight/2 -dy;
	}

	// Drawing line
	ctx.beginPath() ;
	var j ; // creating count variable
	ctx.moveTo(pos_2d[0][0], pos_2d[0][1]);
	for (j=1 ; j<pos_2d.length ; j++){
		ctx.lineTo(pos_2d[j][0], pos_2d[j][1]);
	}
	ctx.stroke();

}

// "right": negative rotation of projection plane around axis OZ
document.getElementById("bt_right").onclick = function() {

	if(document.getElementById("checkbox_translate").checked){
		dx = dx - 5 ; // shifting object from center
	} else {
		beta = beta - 5 * Math.PI/180 // rotating projection plane
	}

	// clearing the canvas
	ctx.clearRect(0, 0, cheight, cheight) ;

	// Projecting points on 2D plane
	var j ; // creating count variable
	for(j=0 ; j<pos.length ; j++){
		pos_2d[j] = (proj(pos[j], alpha, beta, delta)) ;
	}

	// converting positions to fit the HTML canvas frame
	var j ; // creating count variable
	for (j=0 ; j<pos_2d.length ; j++){
		pos_2d[j][0] = pos_2d[j][0]+cheight/2 +dx;
		pos_2d[j][1] = -pos_2d[j][1]+cheight/2 -dy;
	}

	// Drawing line
	ctx.beginPath() ;
	var j ; // creating count variable
	ctx.moveTo(pos_2d[0][0], pos_2d[0][1]);
	for (j=1 ; j<pos_2d.length ; j++){
		ctx.lineTo(pos_2d[j][0], pos_2d[j][1]);
	}
	ctx.stroke();

}

// "clockwise": positive rotation of coordinates along normal to projection plane
document.getElementById("bt_clockwise").onclick = function() {

	delta = delta + 5 * Math.PI/180 // rotating projection plane

	// clearing the canvas
	ctx.clearRect(0, 0, cheight, cheight) ;

	// Projecting points on 2D plane
	var j ; // creating count variable
	for(j=0 ; j<pos.length ; j++){
		pos_2d[j] = (proj(pos[j], alpha, beta, delta)) ;
	}

	// converting positions to fit the HTML canvas frame
	var j ; // creating count variable
	for (j=0 ; j<pos_2d.length ; j++){
		pos_2d[j][0] = pos_2d[j][0]+cheight/2;
		pos_2d[j][1] = -pos_2d[j][1]+cheight/2;
	}

	// Drawing line
	ctx.beginPath() ;
	var j ; // creating count variable
	ctx.moveTo(pos_2d[0][0], pos_2d[0][1]);
	for (j=1 ; j<pos_2d.length ; j++){
		ctx.lineTo(pos_2d[j][0], pos_2d[j][1]);
	}
	ctx.stroke();

}

// "counterclockwise": negative rotation of coordinates along normal to projection plane
document.getElementById("bt_counterclockwise").onclick = function() {

	delta = delta - 5 * Math.PI/180 // rotating projection plane

	// clearing the canvas
	ctx.clearRect(0, 0, cheight, cheight) ;

	// Projecting points on 2D plane
	var j ; // creating count variable
	for(j=0 ; j<pos.length ; j++){
		pos_2d[j] = (proj(pos[j], alpha, beta, delta)) ;
	}

	// converting positions to fit the HTML canvas frame
	var j ; // creating count variable
	for (j=0 ; j<pos_2d.length ; j++){
		pos_2d[j][0] = pos_2d[j][0]+cheight/2;
		pos_2d[j][1] = -pos_2d[j][1]+cheight/2;
	}

	// Drawing line
	ctx.beginPath() ;
	var j ; // creating count variable
	ctx.moveTo(pos_2d[0][0], pos_2d[0][1]);
	for (j=1 ; j<pos_2d.length ; j++){
		ctx.lineTo(pos_2d[j][0], pos_2d[j][1]);
	}
	ctx.stroke();

}

// zoom in
document.getElementById("bt_zoom_in").onclick = function() {

	// increasing values of positions array
	pos = list_prod(pos,1.1) ;

	// clearing the canvas
	ctx.clearRect(0, 0, cheight, cheight) ;

	// Projecting points on 2D plane
	var j ; // creating count variable
	for(j=0 ; j<pos.length ; j++){
		pos_2d[j] = (proj(pos[j], alpha, beta, delta)) ;
	}

	// converting positions to fit the HTML canvas frame
	var j ; // creating count variable
	for (j=0 ; j<pos_2d.length ; j++){
		pos_2d[j][0] = pos_2d[j][0]+cheight/2 +dx;
		pos_2d[j][1] = -pos_2d[j][1]+cheight/2 -dy;
	}

	// Drawing line
	ctx.beginPath() ;
	var j ; // creating count variable
	ctx.moveTo(pos_2d[0][0], pos_2d[0][1]);
	for (j=1 ; j<pos_2d.length ; j++){
		ctx.lineTo(pos_2d[j][0], pos_2d[j][1]);
	}
	ctx.stroke();

}

// zoom out
document.getElementById("bt_zoom_out").onclick = function() {

	// decreasing values of positions array
	pos = list_prod(pos,0.9) ;

	// clearing the canvas
	ctx.clearRect(0, 0, cheight, cheight) ;

	// Projecting points on 2D plane
	var j ; // creating count variable
	for(j=0 ; j<pos.length ; j++){
		pos_2d[j] = (proj(pos[j], alpha, beta, delta)) ;
	}

	// converting positions to fit the HTML canvas frame
	var j ; // creating count variable
	for (j=0 ; j<pos_2d.length ; j++){
		pos_2d[j][0] = pos_2d[j][0]+cheight/2 +dx;
		pos_2d[j][1] = -pos_2d[j][1]+cheight/2 -dy;
	}

	// Drawing line
	ctx.beginPath() ;
	var j ; // creating count variable
	ctx.moveTo(pos_2d[0][0], pos_2d[0][1]);
	for (j=1 ; j<pos_2d.length ; j++){
		ctx.lineTo(pos_2d[j][0], pos_2d[j][1]);
	}
	ctx.stroke();

}

// "adapt size": resizing the drawing to fit the canvas
document.getElementById("bt_adapt_size").onclick = function() {

	// resizing positions to fit the graph
	pos = list_prod(pos, cheight/2/getmax(pos)*9/10) ;

	// clearing the canvas
	ctx.clearRect(0, 0, cheight, cheight) ;

	// Projecting points on 2D plane
	var j ; // creating count variable
	for(j=0 ; j<pos.length ; j++){
		pos_2d[j] = (proj(pos[j], alpha, beta, delta)) ;
	}

	// converting positions to fit the HTML canvas frame
	var j ; // creating count variable
	for (j=0 ; j<pos_2d.length ; j++){
		pos_2d[j][0] = pos_2d[j][0]+cheight/2 ;
		pos_2d[j][1] = -pos_2d[j][1]+cheight/2 ;
	}

	// Drawing line
	ctx.beginPath() ;
	var j ; // creating count variable
	ctx.moveTo(pos_2d[0][0], pos_2d[0][1]);
	for (j=1 ; j<pos_2d.length ; j++){
		ctx.lineTo(pos_2d[j][0], pos_2d[j][1]);
	}
	ctx.stroke();

}

// play animation
document.getElementById("bt_animation").onclick = function() {

	// clearing the canvas
	ctx.clearRect(0, 0, cheight, cheight) ;

	// drawing a new point every millisecond
	function animate(){
		var j = 0 ;
		timer = setInterval(function(){
			ctx.beginPath();
			ctx.moveTo(pos_2d[j][0],pos_2d[j][1]);
			ctx.lineTo(pos_2d[j+1][0],pos_2d[j+1][1]);
			ctx.stroke();
			j++;
		}, 1);
	}
	animate();

}

// stop animation
document.getElementById("bt_animation_stop").onclick = function(){
	clearTimeout(timer);
}




/*Third part of the script*****************************************

This part of the script contains functions that allow the user to
store the drawing, as well as a request to display all images from
the gallery on the page

Functions:
	- Download picture
	- Upload picture to the gallery

*******************************************************************/


// download picture
document.getElementById("bt_dl").onclick = function() {

	/*
	In this function, a new canvas with a different resolution is created
	The trajectory is plotted again on this new canvas, and it is downloaded
	*/

	if(document.getElementById("dl_res").value == 600){

		// creating a 600x600px canvas
		var canvas_tobedl = document.createElement("canvas");
		canvas_tobedl.width = 600;
		canvas_tobedl.height = 600;

	} else if (document.getElementById("dl_res").value == 640){

		// creating a 640x480px canvas
		var canvas_tobedl = document.createElement("canvas");
		canvas_tobedl.width = 640;
		canvas_tobedl.height = 480;

	} else if (document.getElementById("dl_res").value == 1024){

		// creating a 1024x748px canvas
		var canvas_tobedl = document.createElement("canvas");
		canvas_tobedl.width = 1024;
		canvas_tobedl.height = 768;

	} else if (document.getElementById("dl_res").value == 1152){

		// creating a 1152x864px canvas
		var canvas_tobedl = document.createElement("canvas");
		canvas_tobedl.width = 1152;
		canvas_tobedl.height = 864;

	} else if (document.getElementById("dl_res").value == 1600){

		// creating a 1600x1200px canvas
		var canvas_tobedl = document.createElement("canvas");
		canvas_tobedl.width = 1600;
		canvas_tobedl.height = 1200;

	}

	canvas_tobedl_ctx = canvas_tobedl.getContext("2d");
	canvas_tobedl_ctx.lineWidth = 0.3 ;

	// adjusting image size to new canvas resolution
	pos = list_prod(pos, canvas_tobedl.height/cheight) ;

	// Projecting points on 2D plane
	var j ; // creating count variable
	var pos_2d = [[]] ; // creating array of 2D positions
	for(j=0 ; j<pos.length ; j++){
		pos_2d[j] = (proj(pos[j], alpha, beta, delta)) ;
	}

	// converting positions to fit the HTML canvas frame
	var j ; // creating count variable
	for (j=0 ; j<pos_2d.length ; j++){
		pos_2d[j][0] = pos_2d[j][0]+canvas_tobedl.width/2 +dx;
		pos_2d[j][1] = -pos_2d[j][1]+canvas_tobedl.height/2 -dy;
	}

	// Drawing line
	canvas_tobedl_ctx.beginPath() ;
	var j ; // creating count variable
	canvas_tobedl_ctx.moveTo(pos_2d[0][0], pos_2d[0][1]);
	for (j=1 ; j<pos_2d.length ; j++){
		canvas_tobedl_ctx.lineTo(pos_2d[j][0], pos_2d[j][1]);
	}
	canvas_tobedl_ctx.stroke();

	// creating a (fictional) link to download the file
	var download_link = document.createElement("a");
	download_link.download = "Erdorbit.png"; // file name
	download_link.href = canvas_tobedl.toDataURL(); // URL of canvas data
	download_link.click(); // trigger click of the link

}

// upload picture to gallery
document.getElementById("bt_upload").onclick = function() {

	/*
	In this function, the parameters entered by the user are stored
	in a variable, and sent to the server along with the canvas data
	through a php script.
	The php script stores the parameters in a .csv file, and stores
	the canvas data in a .png file in a gallery folder
	*/

	// converting drawing projection plane angles to degrees
	var alpha_deg = alpha * 180/Math.PI ;
	var beta_deg = beta * 180/Math.PI ;
	var delta_deg = delta * 180/Math.PI ;

	// creating string containing drawing parameters
	var drawing_param =
		String(a) + "," +
		String(e) + "," +
		String(i) + "," +
		String(om) + "," +
		String(dur/86400) + "," + // simulation duration is stored in days
		String(step_nb) ;

	// getting username
	var username = document.getElementById("username").value ;
	var useremail = document.getElementById("useremail").value ;

	// getting URL of canvas data
	var canvas_dataURL = canvas.toDataURL();

	// ajax request to send data to server
	$.ajax({
		type: "POST", // type of ajax request best fit for large amount of data
	    url: "upload.php", // URL to send request to
	    data: {
			imgBase64: canvas_dataURL, // sending canvas data encrypted in base 64
			username: username,
			useremail: useremail,
			drawing_param : drawing_param
	  }
	}).done(function() {
	  alert('upload complete'); // displaying success to the user
	  console.log('saved'); // displaying success in the console
	});

}

// displaying all png images from the gallery directory to the html page
$.ajax({ // retrieving database info from server with ajax request
	type: "GET",
	url: "uploads_database.csv", // location of database
	success: function (data) {

		// transferring drawings info from server to javascript array
		var csv = [] ; // array that will receive drawings info
		var data_lines = data.split(/\r\n|\n/); // split lines from csv database
		var headers = data_lines[0].split(','); // extract headers (first line)

		for (var i=1; i<data_lines.length; i++) {
			var data_cell = data_lines[i].split(','); // split cells from csv
			if (data_cell.length == headers.length) {

				var single_line = [];
				for (var j=0; j<headers.length; j++) {
					single_line.push(data_cell[j].slice(1,-1)); // get rid of double quotes
				}
				csv.push(single_line);
			}
		}

		// append drawings info to html file
		for(j=0;j<csv.length;j++){
			if (csv[j][2]){ // in case a username was given
				$('#gallery_div').append("<img src='gallery/"+
				"Bild"+csv[j][0]+"_"+csv[j][1]+"_"+csv[j][4]+".png'><br>"+ // append image link
				"<div class='tooltip'>"+"Bild"+csv[j][0]+"_"+csv[j][1]+ // image label text
				"<span class='tooltiptext'>"+ // tooltip image info
				"<u>Uploaded by:</u><br>"+csv[j][2]+"<br>"+
				"<u>Parameters:</u><br>"+
				"a(km):"+csv[j][5]+"<br>"+
				"e(-):"+csv[j][6]+"<br>"+
				"i(deg):"+csv[j][7]+"<br>"+
				"om(deg):"+csv[j][8]+"<br>"+
				"dur(days):"+csv[j][9]+"<br>"+
				"step(-):"+csv[j][10]+"<br>"+
				"</span></div><br><br>"
				);
			} else { // in cas no username was given
				$('#gallery_div').append("<img src='gallery/"+
				"Bild"+csv[j][0]+"_"+csv[j][1]+"_"+csv[j][4]+".png'><br>"+ // append image link
				"<div class='tooltip'>"+"Bild"+csv[j][0]+"_"+csv[j][1]+ // image label text
				"<span class='tooltiptext'>"+ // tooltip image info
				"<u>Parameters:</u><br>"+
				"a(km):"+csv[j][5]+"<br>"+
				"e(-):"+csv[j][6]+"<br>"+
				"i(deg):"+csv[j][7]+"<br>"+
				"om(deg):"+csv[j][8]+"<br>"+
				"dur(days):"+csv[j][9]+"<br>"+
				"step(-):"+csv[j][10]+"<br>"+
				"</span></div><br><br>"
				);
			}
		}
	}
});
