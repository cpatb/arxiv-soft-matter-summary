"""
ArXiv Soft Matter Paper Summary Generator
Fetches recent papers from cond-mat.soft and generates conceptual synopses
"""

import arxiv
import requests
from datetime import datetime, timedelta
import json
import time
from pathlib import Path

def fetch_recent_papers(days_back=1):
    """Fetch papers from arXiv cond-mat.soft from specified days ago"""

    # Calculate date range
    target_date = datetime.now() - timedelta(days=days_back)
    date_str = target_date.strftime("%Y%m%d")

    # Search query for soft condensed matter
    search = arxiv.Search(
        query="cat:cond-mat.soft",
        max_results=100,
        sort_by=arxiv.SortCriterion.SubmittedDate,
        sort_order=arxiv.SortOrder.Descending
    )

    papers = []
    for result in search.results():
        # Check if paper was submitted on target date
        submission_date = result.published.strftime("%Y%m%d")
        if submission_date == date_str:
            papers.append({
                'title': result.title,
                'authors': [author.name for author in result.authors],
                'abstract': result.abstract,
                'arxiv_id': result.entry_id.split('/')[-1],
                'pdf_url': result.pdf_url,
                'published': result.published.strftime("%Y-%m-%d"),
                'categories': result.categories
            })

    return papers, target_date

def generate_synopsis(title, abstract, api_key=None):
    """Generate a 2-paragraph conceptual synopsis using Perplexity API or fallback"""

    if api_key:
        try:
            headers = {
                "Authorization": f"Bearer {api_key}",
                "Content-Type": "application/json"
            }

            prompt = f"""Write a 2-paragraph conceptual synopsis of this physics paper. 
First paragraph: Explain the main concept, motivation, and research question in accessible terms.
Second paragraph: Summarize the key findings, methodology, and significance.

Title: {title}
Abstract: {abstract}

Synopsis:"""

            payload = {
                "model": "llama-3.1-sonar-small-128k-online",
                "messages": [{"role": "user", "content": prompt}],
                "temperature": 0.3,
                "max_tokens": 500
            }

            response = requests.post(
                "https://api.perplexity.ai/chat/completions",
                json=payload,
                headers=headers,
                timeout=30
            )

            if response.status_code == 200:
                return response.json()['choices'][0]['message']['content']
            else:
                print(f"API error: {response.status_code}")
        except Exception as e:
            print(f"Error generating synopsis: {e}")

    # Fallback: create a basic summary from abstract
    sentences = abstract.replace('\n', ' ').split('. ')
    if len(sentences) > 4:
        para1 = '. '.join(sentences[:len(sentences)//2]) + '.'
        para2 = '. '.join(sentences[len(sentences)//2:]) + '.'
        return f"{para1}\n\n{para2}"
    else:
        return abstract

def generate_html_report(papers, target_date, output_file="arxiv_summary.html"):
    """Generate HTML report with paper summaries"""

    api_key = os.environ.get('PERPLEXITY_API_KEY')

    html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ArXiv Soft Matter Papers - {target_date.strftime('%B %d, %Y')}</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
            line-height: 1.6;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .header h1 {{
            margin: 0 0 10px 0;
            font-size: 2em;
        }}
        .header .date {{
            font-size: 1.2em;
            opacity: 0.9;
        }}
        .paper {{
            background: white;
            padding: 25px;
            margin-bottom: 25px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            transition: transform 0.2s, box-shadow 0.2s;
        }}
        .paper:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.15);
        }}
        .paper-title {{
            color: #2c3e50;
            font-size: 1.4em;
            margin-bottom: 10px;
            font-weight: 600;
        }}
        .paper-title a {{
            color: #667eea;
            text-decoration: none;
        }}
        .paper-title a:hover {{
            text-decoration: underline;
        }}
        .authors {{
            color: #666;
            font-size: 0.95em;
            margin-bottom: 15px;
            font-style: italic;
        }}
        .synopsis {{
            color: #34495e;
            text-align: justify;
            margin-bottom: 15px;
        }}
        .synopsis p {{
            margin-bottom: 12px;
        }}
        .links {{
            margin-top: 15px;
            padding-top: 15px;
            border-top: 1px solid #eee;
        }}
        .links a {{
            display: inline-block;
            margin-right: 15px;
            color: #667eea;
            text-decoration: none;
            font-weight: 500;
        }}
        .links a:hover {{
            text-decoration: underline;
        }}
        .paper-count {{
            background: white;
            padding: 15px;
            border-radius: 8px;
            margin-bottom: 20px;
            text-align: center;
            font-weight: 600;
            color: #2c3e50;
        }}
        .loading {{
            text-align: center;
            color: #666;
            font-style: italic;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>📚 ArXiv Soft Condensed Matter Papers</h1>
        <div class="date">Papers submitted on {target_date.strftime('%B %d, %Y')}</div>
    </div>

    <div class="paper-count">
        Found {len(papers)} papers
    </div>
"""

    for i, paper in enumerate(papers, 1):
        print(f"Processing paper {i}/{len(papers)}: {paper['title'][:60]}...")

        # Generate synopsis
        synopsis = generate_synopsis(paper['title'], paper['abstract'], api_key)

        # Format authors
        authors_str = ", ".join(paper['authors'][:5])
        if len(paper['authors']) > 5:
            authors_str += f" (+{len(paper['authors']) - 5} more)"

        html_template += f"""
    <div class="paper">
        <div class="paper-title">
            <a href="https://arxiv.org/abs/{paper['arxiv_id']}" target="_blank">
                {i}. {paper['title']}
            </a>
        </div>
        <div class="authors">{authors_str}</div>
        <div class="synopsis">
            {synopsis.replace(chr(10) + chr(10), '</p><p>').replace(chr(10), ' ')}
        </div>
        <div class="links">
            <a href="https://arxiv.org/abs/{paper['arxiv_id']}" target="_blank">📄 Abstract</a>
            <a href="{paper['pdf_url']}" target="_blank">📑 PDF</a>
            <a href="https://arxiv.org/abs/{paper['arxiv_id']}" target="_blank">arXiv:{paper['arxiv_id']}</a>
        </div>
    </div>
"""

        # Rate limiting
        if api_key and i < len(papers):
            time.sleep(1)

    html_template += """
</body>
</html>
"""

    # Write to file
    output_path = Path(output_file)
    output_path.write_text(html_template, encoding='utf-8')
    print(f"\nReport generated: {output_path.absolute()}")
    return str(output_path.absolute())

def main():
    print("Fetching recent papers from arXiv cond-mat.soft...")
    papers, target_date = fetch_recent_papers(days_back=1)

    if not papers:
        print(f"No papers found for {target_date.strftime('%Y-%m-%d')}")
        print("This might be expected on weekends or holidays.")
        return

    print(f"Found {len(papers)} papers from {target_date.strftime('%B %d, %Y')}")
    print("Generating synopses and HTML report...\n")

    output_file = generate_html_report(papers, target_date)

    print(f"\n✓ Complete! Open the file in your browser:")
    print(f"  {output_file}")

    # Try to open in browser automatically
    import webbrowser
    try:
        webbrowser.open(f'file://{output_file}')
        print("\n(Opened in default browser)")
    except:
        pass

if __name__ == "__main__":
    main()
