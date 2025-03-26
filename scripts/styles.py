def get_style():
    """Read CSS from external file and wrap in HTML head"""
    with open('scripts/styles.css', 'r') as f:
        css_content = f.read()
    return f"""
<head>
    <title>Sample Name</title>
    <style>
    {css_content}
    </style>
</head>
"""


class HTMLTable:
    def __init__(self, virus_2d_array, header, colorfunc, get_color):
        self.virus_array = virus_2d_array
        self.header = header
        self.pallete = colorfunc
        self.get_color = get_color

    def table_head(self):
        """Generate the table header HTML keeping the exact same output format"""
        th_elements = ''.join(f'<th>{c}</th>' for c in self.header)
        return f'<div class="tableFixHead"><table align="center">\n<thead>\n<tr>{th_elements}</tr></thead>'

    def make_td(self, x):
        color = self.get_color(self.pallete, x) if self.get_color else ''
        attrs = {'style': f'background-color: {color}'} if color else {}
        return Tooltip.create(x, tag="td", **attrs)

    def table_body(self):
        """Generate table body with identical output format"""
        rows = []
        for line in self.virus_array:
            cells = ''.join(self.make_td(i) for i in line)
            rows.append(f'<tr>{cells}</tr>')
        return '<tbody>\n' + '\n'.join(rows) + '</tbody>\n</table></div>'

    def render(self):
        """Maintains the exact same method name and concatenation order"""
        return self.table_head() + self.table_body()


class HTMLDetails:
    def __init__(self, details_dict):
        self.details_dict = details_dict

    def render(self):
        """Produces identical HTML output as the original"""
        return ''.join(
            f'<details><summary>{key}</summary>{content}</details>'
            for key, content in self.details_dict.items()
        )


class HTMLImage:
    def __new__(self, encoded_image, type):
        self.encoded_image = encoded_image
        match type:
            case 'jpg':
                self.type = 'jpg'
            case 'svg':
                self.type = 'image/svg+xml;charset=utf-8'
        self.html_image = f'<div class="image-container"><img alt="" src="data:{self.type};base64,{self.encoded_image}" alt="Ooops! This should have been a picture"/></div> \n'
        return self.html_image


class Header:
    def __new__(cls, text, l=3, tooltip_text=""):
        tag = f'h{l}'
        if tooltip_text:
            return f'<{tag} class="tooltip red-header">{text}<span class="tooltip-text">{tooltip_text}</span></{tag}>'
        return f'<{tag}>{text}</{tag}>'


class Tooltip:
    @staticmethod
    def create(content, tag="span", **attrs):
        """
        Creates HTML elements with tooltips
        :param content: Content with optional CONTAMINATION_INFO
        :param tag: HTML tag to create
        :param attrs: Additional HTML attributes
        """
        if "CONTAMINATION_INFO:" not in str(content):
            attrs_str = ' '.join(f'{k}="{v}"' for k, v in attrs.items())
            return f'<{tag} {attrs_str}>{content}</{tag}>'

        main, tooltip = content.split('CONTAMINATION_INFO:', 1)
        attrs_str = ' '.join(f'{k}="{v}"' for k, v in attrs.items())

        if tag == "td":
            classes = "tooltip tooltip-warning"
            if 'class' in attrs:
                classes = f"{attrs['class']} {classes}"
            style = attrs.get("style", "")
            return f'<td class="{classes}" data-tooltip="{tooltip}" style="{style}">{main}</td>'

        return f'<{tag} class="tooltip" data-tooltip="{tooltip}" {attrs_str}>{main}</{tag}>'
