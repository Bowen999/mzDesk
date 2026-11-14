# mzDesk installer for Windows (PowerShell 5.1 or 7, 64-bit Windows on an AMD64 / Intel CPU).
#
#   irm https://mzdesk.com/install.ps1 | iex
#
# Run it again to update. It needs no administrator rights and writes only to %USERPROFILE%\.mzdesk, plus your
# user PATH so that the `mzdesk` command is found.
# (If you saved this file and PowerShell refuses to run it, start it with:
#  powershell -ExecutionPolicy Bypass -File install.ps1 )
#
# What it does, in order:
#   1. find the system          2. find uv (installs it into .mzdesk\bin if there is none)
#   3. read version.json        4. download the wheel and check its sha256
#   5. make .mzdesk\venv        6. install the wheel and its packages
#   7. write the launcher       8. put the launcher on your PATH
#
# Environment:
#   MZDESK_BASE_URL   where version.json and the wheel come from (default https://mzdesk.com)
#   MZDESK_VERSION    install this version instead of the newest one
#   MZDESK_UPDATING   set to 1 by `mzdesk update`: the last line then says "Updated to <version>"
#
# Not verified on Windows yet.

# Everything runs inside one script block: the script may arrive through `iex`, where `exit` would close the
# user's PowerShell window and top-level variables would stay behind in the session.
& {
    $ErrorActionPreference = 'Stop'          # local to this block
    $ProgressPreference = 'SilentlyContinue' # Windows PowerShell 5.1 downloads are many times slower with the progress bar
    # Windows PowerShell 5.1 may not offer TLS 1.2 by default
    [Net.ServicePointManager]::SecurityProtocol = [Net.ServicePointManager]::SecurityProtocol -bor [Net.SecurityProtocolType]::Tls12

    # Run a program and return its exit code. Windows PowerShell 5.1 turns text a program writes to stderr into errors when
    # ErrorActionPreference is 'Stop' (uv prints its progress there), so it is relaxed for the call and the exit code decides.
    function Invoke-Native {
        param([string]$Exe, [string[]]$Arguments, [switch]$Quiet)
        $saved = $ErrorActionPreference
        $ErrorActionPreference = 'Continue'
        try {
            $global:LASTEXITCODE = 1     # stays 1 if the program cannot be started at all
            if ($Quiet) { & $Exe @Arguments *> $null } else { & $Exe @Arguments | Out-Host }
            return $LASTEXITCODE
        } finally {
            $ErrorActionPreference = $saved
        }
    }

    # The steps are written in mzDesk's blue (Cyan is the closest colour every console has); text and result do not depend on it.
    function Write-Step {
        param([string]$Text)
        Write-Host $Text -ForegroundColor Cyan
    }

    # Download a file with a progress bar in the console (Invoke-WebRequest shows nothing here because the progress bar of Windows PowerShell 5.1
    # makes downloads many times slower). Falls back to Invoke-WebRequest when the stream cannot be read or the output is not a console.
    function Save-Download {
        param([string]$Url, [string]$Path)
        $show = $false
        try { $show = -not [Console]::IsOutputRedirected } catch { $show = $false }
        if ($show) {
            try {
                $request = [Net.HttpWebRequest]::Create($Url)
                $response = $request.GetResponse()
                try {
                    $total = [double]$response.ContentLength
                    $source = $response.GetResponseStream()
                    $target = [IO.File]::Create($Path)
                    try {
                        $buffer = New-Object byte[] 1048576
                        $done = 0.0
                        $last = [DateTime]::UtcNow
                        while (($count = $source.Read($buffer, 0, $buffer.Length)) -gt 0) {
                            $target.Write($buffer, 0, $count)
                            $done += $count
                            if ($total -gt 0 -and (([DateTime]::UtcNow - $last).TotalMilliseconds -gt 200)) {
                                $last = [DateTime]::UtcNow
                                $percent = [int](100 * $done / $total)
                                $filled = [int]($percent / 4)
                                $gauge = ('#' * $filled) + ('-' * (25 - $filled))
                                $mbDone = [math]::Round($done / 1MB)
                                $mbTotal = [math]::Round($total / 1MB)
                                Write-Host ("`r  [$gauge] $percent%  $mbDone of $mbTotal MB") -NoNewline -ForegroundColor Cyan
                            }
                        }
                    } finally {
                        $target.Dispose()
                    }
                } finally {
                    $response.Close()
                }
                Write-Host ''
                return
            } catch {
                Write-Host ''
                if (Test-Path -LiteralPath $Path) { Remove-Item -Force -LiteralPath $Path -ErrorAction SilentlyContinue }
            }
        }
        Invoke-WebRequest -UseBasicParsing -Uri $Url -OutFile $Path
    }

    $step = 'starting'
    $tmp = $null
    try {
        $mz = Join-Path $HOME '.mzdesk'
        $bin = "$mz\bin"
        $venv = "$mz\venv"
        $base = 'https://mzdesk.com'
        if ($env:MZDESK_BASE_URL) { $base = $env:MZDESK_BASE_URL.TrimEnd('/') }

        # ---- 1. system
        $step = 'finding your system'
        if ($env:OS -ne 'Windows_NT') { throw 'this installer is for Windows; on macOS and Linux run: curl -fsSL https://mzdesk.com/install.sh | sh' }
        $cpu = $env:PROCESSOR_ARCHITECTURE
        if ($env:PROCESSOR_ARCHITEW6432) { $cpu = $env:PROCESSOR_ARCHITEW6432 }   # a 32-bit PowerShell on 64-bit Windows
        if ($cpu -ne 'AMD64') { throw "the CPU type $cpu is not supported (64-bit Intel / AMD is)." }
        Write-Host "System: windows $cpu"
        New-Item -ItemType Directory -Force -Path "$mz\cache", $bin | Out-Null
        $tmp = "$mz\cache\install-" + [guid]::NewGuid().ToString('N')
        New-Item -ItemType Directory -Force -Path $tmp | Out-Null

        # ---- 2. uv (it fetches Python and installs the packages)
        $step = 'finding uv'
        $uv = $null
        if (Test-Path -LiteralPath "$bin\uv.exe") {
            $uv = "$bin\uv.exe"
        } elseif (Get-Command uv -ErrorAction SilentlyContinue) {
            $uv = (Get-Command uv).Source
        } else {
            $step = 'installing uv'
            Write-Step 'Installing uv, the tool that sets up Python (astral.sh) ...'
            # UV_UNMANAGED_INSTALL: put uv into this folder and leave the PATH alone. A child PowerShell runs astral's script,
            # so that nothing in it can end this one.
            $env:UV_UNMANAGED_INSTALL = $bin
            try {
                $cmd = '[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; irm https://astral.sh/uv/install.ps1 | iex'
                $code = Invoke-Native 'powershell.exe' @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-Command', $cmd)
                if ($code -ne 0) { throw 'the uv installer failed.' }
            } finally {
                Remove-Item Env:UV_UNMANAGED_INSTALL -ErrorAction SilentlyContinue
            }
            $uv = "$bin\uv.exe"
            if (-not (Test-Path -LiteralPath $uv)) { throw "uv was not installed into $bin." }
        }

        # ---- 3. version.json
        $step = 'reading version.json'
        $info = $null
        try {
            $info = (Invoke-WebRequest -UseBasicParsing -Uri "$base/version.json").Content | ConvertFrom-Json
        } catch {
            throw "cannot download $base/version.json (is the internet reachable, and is the address right?)."
        }
        $version = [string]$info.version
        $wheel = [string]$info.wheel
        $sha = [string]$info.sha256
        if ($env:MZDESK_VERSION -and $env:MZDESK_VERSION -ne $version) {
            # an older (or newer) version than the one version.json describes: its checksum sits next to the wheel
            $version = $env:MZDESK_VERSION
            if ($version -notmatch '^[0-9]+\.[0-9]+\.[0-9]+$') { throw 'MZDESK_VERSION must look like 1.0.0.' }
            $wheel = "mzdesk-$version-py3-none-any.whl"
            try {
                $text = (Invoke-WebRequest -UseBasicParsing -Uri "$base/download/$wheel.sha256").Content
                $sha = ([string]$text).Trim().Split(' ')[0]
            } catch {
                throw "version $version is not available on $base."
            }
        }
        if ($version -notmatch '^[0-9]+\.[0-9]+\.[0-9]+$') { throw 'version.json has no valid version.' }
        if ($wheel -cne "mzdesk-$version-py3-none-any.whl") { throw 'version.json names an unexpected wheel file.' }
        if ($sha -notmatch '^[0-9a-f]{64}$') { throw 'version.json has no valid sha256.' }

        # ---- 4. wheel
        $step = "downloading mzDesk $version"
        Write-Step "Downloading mzDesk $version ..."
        $wheelPath = "$tmp\$wheel"
        try {
            Save-Download -Url "$base/download/$wheel" -Path $wheelPath
        } catch {
            throw "cannot download the wheel from $base/download/$wheel."
        }
        $step = 'checking the download'
        $got = (Get-FileHash -Algorithm SHA256 -LiteralPath $wheelPath).Hash.ToLower()
        if ($got -ne $sha) { throw 'the downloaded file does not match its sha256 checksum; try again, and report it if it repeats.' }

        # ---- 5. Python environment. A working one is kept, so an update only replaces the packages that changed.
        $step = 'creating the Python environment'
        $py = "$venv\Scripts\python.exe"
        $usable = $false
        if (Test-Path -LiteralPath $py) {
            $usable = ((Invoke-Native $py @('-c', 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)') -Quiet) -eq 0)
        }
        if ($usable) {
            Write-Host "Using the Python environment in $venv"
        } else {
            Write-Step "Setting up Python 3.12 in $venv ..."
            if (Test-Path -LiteralPath $venv) { Remove-Item -Recurse -Force -LiteralPath $venv }
            # --managed-python: Python comes from uv, so the environment does not depend on a Python of the system
            $code = Invoke-Native $uv @('venv', $venv, '--python', '3.12', '--managed-python')
            if ($code -ne 0) { throw 'uv could not create the Python environment (it needs the internet to fetch Python).' }
        }

        # ---- 6. packages
        $step = 'installing mzDesk and its packages'
        Write-Step "Installing mzDesk $version and its packages (about 1 GB; a few minutes) ..."
        # --reinstall-package: a rebuilt wheel of the same version replaces the old one; the other packages stay if they still fit
        $code = Invoke-Native $uv @('pip', 'install', '--python', $py, '--reinstall-package', 'mzdesk', $wheelPath)
        if ($code -ne 0) { throw 'uv could not install the packages (see the message above; the internet is needed).' }

        # ---- 7. launcher
        $step = 'writing the launcher'
        $launcher = "@echo off`r`n" +
            "rem Written by the mzDesk installer: runs mzDesk from %USERPROFILE%\.mzdesk\venv and opens it in the browser.`r`n" +
            "`"%USERPROFILE%\.mzdesk\venv\Scripts\python.exe`" -m mzdesk --open %*`r`n"
        Set-Content -LiteralPath "$bin\mzdesk.cmd" -Value $launcher -Encoding ASCII

        # ---- 8. PATH (the user's, no administrator rights)
        $step = 'putting mzdesk on your PATH'
        $newTerminal = $false
        $userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
        $entries = @()
        if ($userPath) { $entries = @($userPath -split ';' | Where-Object { $_ }) }
        if ($entries -notcontains $bin) {
            [Environment]::SetEnvironmentVariable('Path', (($entries + $bin) -join ';'), 'User')
            Write-Host "Added $bin to your PATH"
            $newTerminal = $true
        }
        if (@($env:Path -split ';') -notcontains $bin) { $env:Path = "$env:Path;$bin" }   # this window can use it at once
        if ($newTerminal) { Write-Host 'Terminal windows that are already open must be closed and opened again to find it.' }

        if ($env:MZDESK_UPDATING -eq '1') {
            Write-Step "Updated to $version"
        } else {
            Write-Step "mzDesk $version is installed. Run: mzdesk"
        }
    } catch {
        Write-Host "mzDesk was not installed. Step `"$step`": $($_.Exception.Message)" -ForegroundColor Red
        $global:LASTEXITCODE = 1
    } finally {
        if ($tmp -and (Test-Path -LiteralPath $tmp)) { Remove-Item -Recurse -Force -LiteralPath $tmp -ErrorAction SilentlyContinue }
    }
}
