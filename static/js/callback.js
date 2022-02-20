$('#genes').on('change',function(){
    $.ajax({
        url: "/add_genes",
        type: "GET",
        contentType: 'application/json;charset=UTF-8',
        data: {
            'selected': document.getElementById('genes').value
        },
        dataType:"json",
        success: function (data) {
            Plotly.newPlot('chart', data );
        }
    });
})
