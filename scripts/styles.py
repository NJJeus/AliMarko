import urllib.parse

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

def create_sequence_utility_html(contig_id: str, sequence: str) -> str:
    """
    Generates HTML for a sequence viewer, copy-to-clipboard button, and download button.
    """
    if not sequence:
        return f'<p class="text-red-500">Error: No sequence found for {contig_id}</p>'

    # URL-encode the sequence string for the data URI in the download link
    encoded_sequence = urllib.parse.quote(f'>{contig_id}\n{sequence}')

    # HTML structure with Tailwind classes (assuming the styles include Tailwind)
    html_content = f"""
    <div class="p-4 mb-4 bg-gray-50 rounded-lg shadow-inner border border-gray-200">
        <h4 class="text-lg font-semibold text-gray-700 mb-2">Sequence ({len(sequence):,} bp)</h4>
        
        <!-- Hidden Textarea to hold the FASTA sequence for reliable copying -->
        <textarea id="seq-content-{contig_id}" class="hidden" style="visibility: hidden; height: 0px;">>{contig_id}\n{sequence}</textarea>

        
        <!-- Buttons -->
        <div class="flex flex-col sm:flex-row gap-3">
            
            <!-- Download Button (Pure HTML Anchor Tag with data URI) -->
            <a href="data:text/fasta;charset=utf-8,{encoded_sequence}" 
               download="{contig_id}.fasta" 
               class="flex-1 text-center px-4 py-2 text-sm font-medium rounded-lg shadow-md text-white bg-indigo-600 hover:bg-indigo-700 transition duration-150 ease-in-out cursor-pointer">
                Download FASTA ({contig_id}.fasta)
            </a>

            <!-- Copy to Clipboard Button (Requires minimal JS) -->
            <button onclick="copySequenceToClipboard('{contig_id}')" 
                    class="flex-1 text-center px-4 py-2 text-sm font-medium rounded-lg shadow-md text-white bg-green-600 hover:bg-green-700 transition duration-150 ease-in-out">
                Copy to Clipboard
            </button>
        </div>
    </div>
    """
    return html_content



def get_utility_script() -> str:
    """
    Returns the JavaScript needed for the copy-to-clipboard functionality,
    including a custom message box instead of alert().
    """
    # Note: This is embedded at the end of the <body> in generate_html_output
    # It uses a global message box defined in the CSS (which is usually part of styles.get_style())
    # I am embedding the message box element and the required JS here for a complete solution.
    return """
<div id="message-box" class="fixed top-4 left-1/2 transform -translate-x-1/2 p-3 rounded-lg text-white bg-gray-800 shadow-xl opacity-0 transition-opacity duration-300 z-50 pointer-events-none">
    Copied successfully!
</div>
<style>
    /* Ensure the message box is hidden initially and positioned correctly */
    #message-box {
        opacity: 0;
        visibility: hidden;
    }
    #message-box.visible {
        opacity: 1;
        visibility: visible;
    }
</style>
<script>
    /**
     * Shows a transient success message box instead of using alert().
     * @param {string} message - The message to display.
     */
    function showMessage(message) {
        const msgBox = document.getElementById('message-box');
        msgBox.textContent = message;
        msgBox.classList.add('visible');

        setTimeout(() => {
            msgBox.classList.remove('visible');
        }, 2000); // Display for 2 seconds
    }

    /**
     * Copies the content of the hidden textarea corresponding to the contig ID
     * to the clipboard.
     * @param {string} contigId - The ID of the contig (used to find the textarea).
     */
    function copySequenceToClipboard(contigId) {
        const textarea = document.getElementById(`seq-content-${contigId}`);
        if (!textarea) {
            showMessage('Copy failed: Sequence content not found.');
            return;
        }

        const content = textarea.value;

        // Use the modern Clipboard API
        if (navigator.clipboard && navigator.clipboard.writeText) {
            navigator.clipboard.writeText(content)
                .then(() => {
                    showMessage('FASTA copied to clipboard!');
                })
                .catch(err => {
                    console.error('Could not copy text: ', err);
                    showMessage('Failed to copy. Check console for details.');
                });
        } else {
            // Fallback using document.execCommand (less reliable in modern browsers/environments)
            try {
                textarea.select();
                document.execCommand('copy'); 
                showMessage('FASTA copied to clipboard via fallback!');
            } catch (err) {
                console.error('Fallback copy failed: ', err);
                showMessage('Copy failed (no clipboard support).');
            }
        }
    }
</script>
"""

