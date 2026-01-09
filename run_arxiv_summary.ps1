# ArXiv Daily Summary - PowerShell Launcher
# Double-click this file or run from PowerShell to get your daily papers

Write-Host "ArXiv Soft Matter Daily Summary" -ForegroundColor Cyan
Write-Host "================================" -ForegroundColor Cyan
Write-Host ""

# Check if Python is available
try {
    $pythonCmd = Get-Command python -ErrorAction Stop
    Write-Host "✓ Python found: $($pythonCmd.Source)" -ForegroundColor Green
} catch {
    Write-Host "✗ Python not found in PATH" -ForegroundColor Red
    Write-Host "Please install Python or ensure it's in your PATH" -ForegroundColor Yellow
    Read-Host "Press Enter to exit"
    exit 1
}

# Run the script
Write-Host ""
Write-Host "Fetching papers and generating report..." -ForegroundColor Yellow
Write-Host ""

python arxiv_daily_summary.py

Write-Host ""
Write-Host "Press Enter to close..." -ForegroundColor Gray
Read-Host
