var $map = jQuery.noConflict();

/* bodymap*/

$map(".pic_1").hover(function(e) {
  var id=$map(this).attr("id");
  //alert(id);
  var strArray = id.split("_");
  //alert(strArray[0]);
  $map("#"+strArray[0]+"_pic").addClass("hover-item");
  $map("#"+strArray[0]+"_text").addClass("select_text");
  
  $map("#"+strArray[0]+"_pic_table").addClass("hover_table_data");
  //alert(id); 
});

$map(".pic_1").mouseleave(function(e) {
  var id=$map(this).attr("id");
  //alert(id);
  var strArray = id.split("_");
  $map("#"+strArray[0]+"_pic").removeClass("hover-item");
  $map("#"+strArray[0]+"_text").removeClass("select_text");
  $map("#"+strArray[0]+"_pic_table").removeClass("hover_table_data");
});

$map(".pic").hover(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  //alert(strArray[0]);
	  $map("#"+strArray[0]+"_node").addClass("hover-item");
	  $map("#"+strArray[0]+"_text").addClass("select_text");
	  $map("#"+strArray[0]+"_pic_table").addClass("hover_table_data");
	  //alert("#"+strArray[0]+"_name");
});

$map(".pic").mouseleave(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  //alert(strArray[0]);
	  $map("#"+strArray[0]+"_node").removeClass("hover-item");
	  $map("#"+strArray[0]+"_text").removeClass("select_text");
	  $map("#"+strArray[0]+"_pic_table").removeClass("hover_table_data");
	  //alert(id);
});


$map(".pic_table").hover(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  //alert(strArray[0]);
	  $map("#"+strArray[0]+"_node").addClass("hover-item");
	  $map("#"+strArray[0]+"_text").addClass("select_text");
	  $map("#"+strArray[0]+"_pic").addClass("hover-item");
	  //alert("#"+strArray[0]+"_name");
});

$map(".pic_table").mouseleave(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  //alert(strArray[0]);
	  $map("#"+strArray[0]+"_node").removeClass("hover-item");
	  $map("#"+strArray[0]+"_text").removeClass("select_text");
	  $map("#"+strArray[0]+"_pic").removeClass("hover-item");
	  //alert(id);
});




$map(".pic_text").hover(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  //alert(strArray[0]);
	  $map("#"+strArray[0]+"_node").addClass("hover-item");
	  $map("#"+strArray[0]+"_pic").addClass("hover-item");
	  $map("#"+strArray[0]+"_pic_table").addClass("hover_table_data");
	  //alert("#"+strArray[0]+"_name");
});

$map(".pic_text").mouseleave(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  //alert(strArray[0]);
	  $map("#"+strArray[0]+"_node").removeClass("hover-item");
	  $map("#"+strArray[0]+"_pic").removeClass("hover-item");
	  $map("#"+strArray[0]+"_pic_table").removeClass("hover_table_data");
	  //alert(id);
});

/*  cellmap  1*/


$map(".cell_dark").hover(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
//	  alert(strArray[0]+"_pic");
	  $map("#"+strArray[0]+"_pic").addClass("hover_cell");
	  $map("#"+strArray[0]+"_text").addClass("hover_text");
	  $map("#"+strArray[0]+"_icon").addClass("hover_icon");
	  //alert(id); 
	});

$map(".cell_dark").mouseleave(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  $map("#"+strArray[0]+"_pic").removeClass("hover_cell");
	  $map("#"+strArray[0]+"_text").removeClass("hover_text");
	  $map("#"+strArray[0]+"_icon").removeClass("hover_icon");
	});


$map(".cell_text").hover(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
//	  alert(strArray[0]+"_pic");
	  $map("#"+strArray[0]+"_pic").addClass("hover_cell");
	  $map("#"+strArray[0]+"_text").addClass("hover_text");
	  $map("#"+strArray[0]+"_icon").addClass("hover_icon");
	  //alert(id); 
	});

$map(".cell_text").mouseleave(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  $map("#"+strArray[0]+"_pic").removeClass("hover_cell");
	  $map("#"+strArray[0]+"_text").removeClass("hover_text");
	  $map("#"+strArray[0]+"_icon").removeClass("hover_icon");
	});


$map(".cell_icon").hover(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
//	  alert(strArray[0]+"_pic");
	  $map("#"+strArray[0]+"_pic").addClass("hover_cell");
	  $map("#"+strArray[0]+"_text").addClass("hover_text");
	  $map("#"+strArray[0]+"_icon").addClass("hover_icon");
	  //alert(id); 
	});

$map(".cell_icon").mouseleave(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  $map("#"+strArray[0]+"_pic").removeClass("hover_cell");
	  $map("#"+strArray[0]+"_text").removeClass("hover_text");
	  $map("#"+strArray[0]+"_icon").removeClass("hover_icon");
	});



/* Muti text 细胞核*/

$map(".sub_cell_dark").hover(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
//	  alert(strArray[0]+"_pic");
	  $map("#"+strArray[0]+"_pic").addClass("hover_cell");
	  $map("#"+strArray[0]+"_text").addClass("subtext");
	  $map("#Nucleolus").addClass("outborder");
	  $map("#Nucleolus_icon").addClass("hover_icon");
	  $map("#all_title").addClass("subtext");
	});

$map(".sub_cell_dark").mouseleave(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  $map("#"+strArray[0]+"_pic").removeClass("hover_cell");
	  $map("#"+strArray[0]+"_text").removeClass("subtext");
	  $map("#Nucleolus").removeClass("outborder");
	  $map("#Nucleolus_icon").removeClass("hover_icon");
	  $map("#all_title").removeClass("subtext");
	});



$map(".small_text").hover(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
//	  alert(strArray[0]+"_pic");
	  $map("#"+strArray[0]+"_pic").addClass("hover_cell");
	  $map("#"+strArray[0]+"_text").addClass("subtext");
	  $map("#Nucleolus").addClass("outborder");
	  $map("#Nucleolus_icon").addClass("hover_icon");
	  $map("#all_title").addClass("subtext");
	});

$map(".small_text").mouseleave(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  $map("#"+strArray[0]+"_pic").removeClass("hover_cell");
	  $map("#"+strArray[0]+"_text").removeClass("subtext");
	  $map("#Nucleolus").removeClass("outborder");
	  $map("#Nucleolus_icon").removeClass("hover_icon");
	  $map("#all_title").removeClass("subtext");
	});



$map(".all_icon").hover(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
//	  alert(strArray[0]+"_pic");
	  $map("#all_title").addClass("subtext");
	  $map("#Nucleolus").addClass("outborder");
	  $map("#Nucleolus_icon").addClass("hover_icon");
	});

$map(".all_icon").mouseleave(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  $map("#all_title").removeClass("subtext");
	  $map("#Nucleolus").removeClass("outborder");
	  $map("#Nucleolus_icon").removeClass("hover_icon");
	});


/* Muti text 外泌体 */


$map(".sub_cell_dark_2").hover(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
//	  alert(strArray[0]+"_pic");
	  $map("#"+strArray[0]+"_pic").addClass("hover_cell");
	  $map("#"+strArray[0]+"_text").addClass("subtext");
	  $map("#Extracellular").addClass("outborder");
	  $map("#EX_icon").addClass("hover_icon");
	  $map("#all_title_2").addClass("subtext");
	});

$map(".sub_cell_dark_2").mouseleave(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  $map("#"+strArray[0]+"_pic").removeClass("hover_cell");
	  $map("#"+strArray[0]+"_text").removeClass("subtext");
	  $map("#Extracellular").removeClass("outborder");
	  $map("#EX_icon").removeClass("hover_icon");
	  $map("#all_title_2").removeClass("subtext");
	});



$map(".small_text_2").hover(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
//	  alert(strArray[0]+"_pic");
	  $map("#"+strArray[0]+"_pic").addClass("hover_cell");
	  $map("#"+strArray[0]+"_text").addClass("subtext");
	  $map("#Extracellular").addClass("outborder");
	  $map("#EX_icon").addClass("hover_icon");
	  $map("#all_title_2").addClass("subtext");
	});

$map(".small_text_2").mouseleave(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  $map("#"+strArray[0]+"_pic").removeClass("hover_cell");
	  $map("#"+strArray[0]+"_text").removeClass("subtext");
	  $map("#Extracellular").removeClass("outborder");
	  $map("#EX_icon").removeClass("hover_icon");
	  $map("#all_title_2").removeClass("subtext");
	});



$map(".all_icon_2").hover(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
//	  alert(strArray[0]+"_pic");
	  $map("#all_title_2").addClass("subtext");
	  $map("#Extracellular").addClass("outborder");
	  $map("#EX_icon").addClass("hover_icon");
	});

$map(".all_icon_2").mouseleave(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  $map("#all_title_2").removeClass("subtext");
	  $map("#Extracellular").removeClass("outborder");
	  $map("#EX_icon").removeClass("hover_icon");
	});



/*  cellmap  2*/


$map(".areatxt").hover(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
//	  alert(strArray[0]+"_pic");
	  $map("#"+strArray[0]+"_pic").addClass("hover_area");
//	  $map("#"+strArray[0]+"_text").addClass("select_text");
	  //alert(id); 
	});

$map(".areatxt").mouseleave(function(e) {
	  var id=$map(this).attr("id");
	  //alert(id);
	  var strArray = id.split("_");
	  $map("#"+strArray[0]+"_pic").removeClass("hover_area");
//	  $map("#"+strArray[0]+"_text").removeClass("select_text");
	});

//	$map(".pic").hover(function(e) {
//		  var id=$map(this).attr("id");
//		  //alert(id);
//		  var strArray = id.split("_");
//		  //alert(strArray[0]);
//		  $map("#"+strArray[0]+"_node").addClass("hover-item");
//		  $map("#"+strArray[0]+"_text").addClass("select_text");
//		  //alert("#"+strArray[0]+"_name");
//	});
//
//	$map(".pic").mouseleave(function(e) {
//		  var id=$map(this).attr("id");
//		  //alert(id);
//		  var strArray = id.split("_");
//		  //alert(strArray[0]);
//		  $map("#"+strArray[0]+"_node").removeClass("hover-item");
//		  $map("#"+strArray[0]+"_text").removeClass("select_text");
//		  //alert(id);
//	});
