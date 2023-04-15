var $map = jQuery.noConflict();
$map(".show").hover(function(e) {
  $map('#info-box').css('display','block');
  $map('#info-box').html('<div></div>'+$map(this).data('info'));
  var  id=$map(this).attr("id");
  //alert(id);
  var classes='.'+id;
  $map(classes).addClass("hover-row");
});

$map(".show").mouseleave(function(e) {
  $map('#info-box').css('display','none');
  var  id=$map(this).attr("id");
  //alert(id);
  var classes='.'+id;
  $map(classes).removeClass("hover-row");
});

$map("#humanbody").mousemove(function(e) {
  //$map('#info-box').css('top',e.pageY-$map('#info-box').height()-210);
  //$map('#info-box').css('left',e.pageX-($map('#info-box').width())+28);
  var isFirefox=navigator.userAgent.toUpperCase().indexOf("CHROME");
  if(isFirefox==-1){
		e = e || window.event;
		//x=e.offsetX;
	    //y=e.offsetY;
		x2=e.clientX;
		y2=e.clientY;
		//var x=$map("#humanbody").offset().left;
		//var y=$map("#humanbody").offset().top;
		//alert(x);
		//alert(y);
	}else{
		//x=event.offsetX;
		//y=event.offsetY;
		x2=event.clientX;
		y2=event.clientY;
		//var x=$map("#humanbody").offset().left;
		//var y=$map("#humanbody").offset().top;
		//alert(x);
		//alert(y);
	}
  $map('#info-box').css('top',e.pageY - $map("div#humanbody").offset().top-$map('#info-box').height()-10);
  $map('#info-box').css('left',e.pageX - $map("div#humanbody").offset().left);
}).mouseover();

$map(".show").click(function() {
//alert($(this).data('disease'));
//window.open("browse_bodymap_res.jsp?type=tissue&searchname="+$map(this).data('disease'),"_self");    ////////////////////my///////////////////////
});


