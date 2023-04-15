///////////////////////////////////////////////
/// 树型标签
/// V 1.0
/// creat by lee
/// https://github.com/miracleren/tagTree
/// 20190529
/// 运行库 juqery
/// //////////////////////////////////////////


;(function($){

	var defaults ={
		id:"",
		data:[],
		fold:true,
		multiple:false,
		check:function(){},
		done:function(){}
	};

	$.fn.tagTree = function(options){
		var that = $(this);
		options.id = "#" + that.attr("id");
        var opts = $.extend(defaults, options);

        that.addClass("tagtree");
        setTree(defaults.data,that);

        $(defaults.id+' li:has(ul)').addClass('li-top');
        if(defaults.fold)
        	$(defaults.id+" .li-top li").hide('fast');
	    $(defaults.id+' li.li-top > span').on('click', function (e) {
	        var children = $(this).parent('li.li-top').find(' > ul > li');
	        if (children.is(":visible")) {
	            children.hide('fast');
	        } else {
	            children.show('fast');
	        }
	        return false;
	    });
	    $(defaults.id+' li span').hover(function() {
	    	if (!$(this).find('i').hasClass('i-check'))
	    		$(this).find('i').show(200);
	    }, function() {
	    	if (!$(this).find('i').hasClass('i-check'))
	    		$(this).find('i').hide(100);
	    });
	    $(defaults.id+' li span i').click(function(event) {
	    	if(!defaults.multiple)
	    	{
	    		$(defaults.id +" .i-check").removeClass('i-check').hide(100);
	    	}

	    	if($(this).hasClass('i-check'))
	    		$(this).removeClass('i-check');
	    	else
	    		$(this).addClass('i-check');
	    	defaults.check($(this).attr("data-val"));
	    	return false;
	    });

	    defaults.done();
	}

	$.fn.tagTreeValues =function(){
		var vals = [];
		$(defaults.id +" .i-check").each(function(index, el) {
			vals.push($(el).attr('data-val'));
		});
		return vals;
	}

	//递归生成树
	function setTree(data,that)
	{
		var ul = $('<ul></ul>');
		that.append(ul);
		$.each(data,function(index,value){
			var str = value.id;
			
				if(str==("GSE_iframe")){
					var li = $('<li><span><a target="right_iframe" style="display:block" href="exmdb_browse_right_iframe_GSE.jsp?searchname='+value.value+'">'+value.name+'<i data-val="'+value.value+'" class="fa fa-check" style="display:none;"></i></span></a></li>');
                }else if(str==("exp_dis")){
                	var li = $('<li><span><a target="right_iframe" style="display:block;text-transform: capitalize;" href="exmdb_exp_data_browse_table.jsp?searchname='+value.value+'">'+value.name+'<i data-val="'+value.value+'" class="fa fa-check" style="display:none;"></i></span></a></li>');
                }else if(str==("marker")){
                	var li = $('<li><span><a target="right_iframe" style="display:block" href="exmdb_browse_biomarker_data_table.jsp?searchname='+value.value+'">'+value.name+'<i data-val="'+value.value+'" class="fa fa-check" style="display:none;"></i></span></a></li>');
                }else if(str==("sc_data")){
                	var li = $('<li><span><a target="right_iframe" style="display:block;" href="exmdb_browse_right_iframe_SC.jsp?searchname='+value.value+'">'+value.name+'<i data-val="'+value.value+'" class="fa fa-check" style="display:none;"></i></span></a></li>');
                }else{
                	var li = $('<li><span><a target="right_iframe" id="qi_va" style="display:block;font-weight: 700;font-size: large;" href="exmdb_exp_data_table.jsp?searchname='+value.value+'">'+value.name+'<i data-val="'+value.value+'" class="fa fa-check" style="display:none;"></i></span></a></li>');
                }
			
				ul.append(li);
		    if(value.children.length > 0)
		    {
		    	li.append('<div class="node-count">'+value.children.length+'</div>');
		    	setTree(value.children,li);
		    }
		});
	}
})(jQuery);