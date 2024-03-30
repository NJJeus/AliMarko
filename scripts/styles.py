

set_style = """
<head>
	<title>Sample Name</title>
	<style>
        .tableFixHead {
        overflow-y: auto; /* make the table scrollable if height is more than 200 px  */
        max-height: 800px; /* gives an initial height of 200px to the table */
        width: 80%;
        margin:auto;
        text-align: center;
                      }
                      
        
        .tableFixHead thead th {
        position: sticky; /* make the table heads sticky */
        top: 0px; /* table head will be placed from the top of the table and sticks to it */
        text-align: center;
        fontsize: 8
                                }
                                
                                
        body {
			background-color: #f1f1f1;
			font-family: Arial, sans-serif;
			color: #333;
			margin: 0;
			padding: 0;
            text-align: center;
		}
		header {
			background-color: #2ecc71;
			color: #000;
			padding: 240px;
            text-align: center;
		}
		h1 {
			font-size: 36px;
			margin-top: 0;
            color: #454444;
            text-align: center;
		}
        h2 {
			font-size: 28px;
			margin-top: 0;
            color: #454444;
            text-align: center;
        }
        h3 {
        font-size: 18px;
        margin-top: 0;
        color: #454444;
        text-align: center;
		}
		main {
			max-width: 1000px;
			margin: 0 auto;
			padding: 20px;
            text-align: center;
		}
		p {
			line-height: 1.5;
			margin-bottom: 20px;
		}
		a {
			color: #2ecc71;
			text-decoration: none;
		}
		a:hover {
			text-decoration: underline;
		}
		table {
			border-collapse: collapse;
			width: 100%;
			margin-bottom: 20px;
			background-color: #fff;
            font-size: 12px;
			box-shadow: 0px 0px 8px rgba(0, 0, 0, 0.1);
            text-align: center;
		}
		th, td {
			text-align: center;
			padding: 5px;
            height: 8px;
            font-size: 12px;
			border-bottom: 1px solid #ddd;
		}
		th {
			background-color: #18aa9d;
			color: #fff;
			font-weight: bold;
            font-size: 12px;
            
            
			text-transform: uppercase;
		}
		tr:hover {
			background-color: #f5f5f5;
		}
        .container {
        position: relative;
        display: inline-block;
        width: 60%;
        min-width: 700px;
    }
    .open-button {
    position: absolute;
    top: 3;
    left: 4;
    padding: 5px;
    background: white;
    cursor: pointer;
    border: 1px solid #ccc;
    border-radius: 5px;
}
	</style>
</head>
"""


class Table:
    def __init__(self, virus_2d_array, header, colorfunc, get_color):
        self.virus_array = virus_2d_array
        self.header = header
        self.pallete = colorfunc
        self.get_color = get_color
    def table_head(self):
        table_head = f'<div class="tableFixHead"><table align="center"> \n <thead> \n <tr>'

        
        table_head += '\n'.join([f'<th>{c}</th>' for c in self.header])
        table_head + '</tr></thead>'

        return table_head

    
    def table_body(self):

        body = '<tbody> \n '
        for line in self.virus_array:
            table_line = '<tr>' + ''.join([f'<td style="background-color: {self.get_color(self.pallete, i)}; width=10; heigh=5">{i}</td>' for i in line]) + '</tr>'
            body += table_line
        body += '</tbody> \n </table></div>'
        return body
    def make_table(self):
        return self.table_head() + self.table_body()
    
class Details:
    def __init__(self, details_dict):
        self.details_dict = details_dict
    def make_details(self):
        details = ''
        for key, content in self.details_dict.items():
            details += f'<details><summary>{key}</summary>'
            details += str(content)
            details += '</details>'
        return details
    
        
        
    

            
        
        
    


