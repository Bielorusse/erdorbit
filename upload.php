<?php

/*
upload.php for Erdorbit
author: thibaut voirand
*/


/*

This script stores a drawing's data in a png file in a gallery folder
on the server, and stores the drawing's parameters in a csv file

This script keeps count of all the images in the gallery folder

*/

// creating alternative to 'file_put_contents' function
// which writes contents in a file
if (!function_exists('file_put_contents')) {
    function file_put_contents($filename, $data) {
        $f = @fopen($filename, 'w'); //@ to hide error, w is the method
        if (!$f) {
            return false;
        } else {
            $bytes = fwrite($f, $data);
            fclose($f);
            return $bytes;
        }
    }
}

// creating alternative to 'fputcsv' function
if (!function_exists('fputcsv')){
	function fputcsv($fh, $arr)
	 {
			$csv = "";
			while (list($key, $val) = each($arr))
			 {
					$val = str_replace('"', '""', $val);
					$csv .= '"'.$val.'",';
			 }
			$csv = substr($csv, 0, -1);
			$csv .= "\n";
			if (!@fwrite($fh, $csv))
			 return FALSE;
	 }
}

// uploading image data, which is crypted in base 64, to the gallery directory

define('UPLOAD_DIR', 'gallery/'); // upload directory

// extracting image data
$img = $_POST['imgBase64'];
$img = str_replace('data:image/png;base64,', '', $img);
$img = str_replace(' ', '+', $img);
$data = base64_decode($img);

// counting images in gallery
$filecount = 0 ;
$files = glob("gallery/" . "*.png");
if ($files){
	$filecount = count($files) + 1 ;
} else {
	$filecount = 1 ; // to start the count at 1 (avoid counting 0,2,3,..)
}

// creating the file name (path, count, date, id, extension)
$random_id = uniqid('',false);
$date = date("Y.m.d");
$file = UPLOAD_DIR . 'Bild' . $filecount . '_' . $date . '_' . $random_id . '.png';
// writing data in the file
$success = file_put_contents($file, $data);
print $success ? $file : 'Unable to save the file.';

// writing upload info in the csv file (database)
$database = fopen("uploads_database.csv","a");
$list = array_merge(array($filecount,$date,$_POST["username"],$_POST["useremail"],$random_id),explode(",",$_POST["drawing_param"]));
fputcsv($database,$list);
fclose($database);

?>
