# ArXiv Soft Matter Daily Summary 📚

Automatically fetch and summarize recent papers from arXiv's Soft Condensed Matter category (cond-mat.soft). Perfect for staying up-to-date with the latest research in your morning routine.

## Features

- 🔍 **Automatic Paper Fetching**: Retrieves papers from the previous day
- 📝 **AI-Generated Synopses**: 2-paragraph summaries explaining motivation, methods, and findings
- 🎨 **Beautiful HTML Reports**: Clean, readable format with gradient headers
- 🔗 **Direct Links**: Quick access to abstracts and PDFs
- ⚡ **Fast Processing**: Typically completes in 1-2 minutes
- 🖱️ **One-Click Execution**: Double-click batch file to run

## Quick Start

### Installation

1. **Clone this repository:**
   ```bash
   git clone https://github.com/yourusername/arxiv-soft-matter-summary.git
   cd arxiv-soft-matter-summary
   ```

2. **Install dependencies:**
   ```bash
   pip install arxiv requests
   ```

### Usage

**Windows - Double-click:**
- Double-click `run_arxiv_summary.bat` or `run_arxiv_summary.ps1`

**Command Line:**
```bash
python arxiv_daily_summary.py
```

The script will automatically:
1. Fetch papers submitted yesterday
2. Generate conceptual synopses
3. Create `arxiv_summary.html`
4. Open the report in your browser

## Example Output

Each paper gets a formatted entry with:
- **Title** (linked to arXiv)
- **Authors**
- **2-Paragraph Synopsis** explaining the research conceptually
- **Links** to abstract and PDF

![Example Screenshot](screenshot.png)

## Optional: Enhanced Synopses

For higher-quality AI-generated synopses, use the Perplexity API:

1. Get an API key from [Perplexity](https://www.perplexity.ai/settings/api)
2. Set environment variable:

**Windows PowerShell:**
```powershell
$env:PERPLEXITY_API_KEY="your-api-key-here"
```

**Linux/Mac:**
```bash
export PERPLEXITY_API_KEY="your-api-key-here"
```

**Permanent (Windows):**
```powershell
[System.Environment]::SetEnvironmentVariable('PERPLEXITY_API_KEY', 'your-api-key-here', 'User')
```

Without the API key, the script generates summaries from abstracts (still useful!).

## Automation

### Windows Task Scheduler

1. Open Task Scheduler → Create Basic Task
2. Name: "ArXiv Daily Summary"
3. Trigger: Daily at 8:00 AM
4. Action: Start a program
   - Program: `python`
   - Arguments: `C:\path\to\arxiv_daily_summary.py`
   - Start in: `C:\path\to\arxiv-soft-matter-summary\`

### Linux/Mac Cron

Add to crontab (`crontab -e`):
```bash
0 8 * * * cd /path/to/arxiv-soft-matter-summary && python3 arxiv_daily_summary.py
```

### Desktop Shortcut

Create a shortcut to `run_arxiv_summary.bat` on your desktop for quick one-click access.

## Customization

Edit `arxiv_daily_summary.py` to customize:

```python
# Fetch papers from N days ago
papers, target_date = fetch_recent_papers(days_back=1)

# Change output filename
output_file = "my_custom_name.html"

# Modify HTML styling
# Edit the CSS in generate_html_report()
```

## Troubleshooting

**"No papers found"**
- ArXiv doesn't receive submissions every day (weekends, holidays)
- Try running later in the day if papers haven't posted yet

**"Python not found"**
- Ensure Python is installed and in your PATH
- Try `python3` instead of `python` on Linux/Mac

**Import errors**
- Run: `pip install arxiv requests`

**Rate limiting**
- The script includes 1-second delays between API calls
- If using Perplexity API, check your rate limits

## File Structure

```
arxiv-soft-matter-summary/
├── arxiv_daily_summary.py    # Main script
├── run_arxiv_summary.bat     # Windows batch launcher
├── run_arxiv_summary.ps1     # PowerShell launcher
├── README.md                 # This file
├── .gitignore               # Git ignore rules
└── arxiv_summary.html       # Generated output (not tracked)
```

## Requirements

- Python 3.7+
- `arxiv` package
- `requests` package

## Contributing

Contributions welcome! Feel free to:
- Report bugs
- Suggest features
- Submit pull requests

## License

MIT License - feel free to use and modify!

## Acknowledgments

- [arXiv](https://arxiv.org/) for providing free access to research papers
- [arxiv-py](https://github.com/lukasschwab/arxiv.py) for the Python API wrapper

---

**Happy reading!** 🔬📖
