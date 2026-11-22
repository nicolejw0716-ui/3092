import { WebR } from 'webr';
import './style.css';

const R_CODE = String.raw`
# Read the uploaded data and keep the original row count for drop reporting.
raw <- read.csv('/input.csv', stringsAsFactors = FALSE, check.names = FALSE)
required <- c('condition', 'substrate_mM', 'velocity')
missing <- setdiff(required, names(raw))
if (length(missing) > 0) stop(paste('Missing required column(s):', paste(missing, collapse = ', ')))

condition_raw <- trimws(as.character(raw$condition))
substrate <- suppressWarnings(as.numeric(as.character(raw$substrate_mM)))
velocity <- suppressWarnings(as.numeric(as.character(raw$velocity)))
valid <- !is.na(condition_raw) & nzchar(condition_raw) & is.finite(substrate) &
  is.finite(velocity) & substrate > 0

dropped_total <- sum(!valid)
d <- data.frame(condition = condition_raw[valid], substrate_mM = substrate[valid], velocity = velocity[valid])
conditions <- unique(condition_raw[!is.na(condition_raw) & nzchar(condition_raw)])
read_counts <- data.frame(
  condition = conditions,
  rows_read = vapply(conditions, function(x) sum(!is.na(condition_raw) & condition_raw == x), integer(1)),
  rows_used = vapply(conditions, function(x) sum(d$condition == x), integer(1))
)
write.csv(read_counts, '/counts.csv', row.names = FALSE, na = '')

# Fit v = Vmax*S/(Km+S) independently for every condition. No transformed fits or pooling.
fits <- lapply(conditions, function(cond) {
  one <- d[d$condition == cond, , drop = FALSE]
  if (nrow(one) < 3) return(list(condition = cond, n = nrow(one), model = NULL, error = 'Need at least 3 valid rows to estimate two parameters.'))
  tryCatch({
    model <- nls(velocity ~ Vmax * substrate_mM / (Km + substrate_mM), data = one,
      start = list(Vmax = max(one$velocity), Km = median(one$substrate_mM)),
      control = nls.control(maxiter = 200, warnOnly = TRUE))
    list(condition = cond, n = nrow(one), model = model, error = '')
  }, error = function(e) list(condition = cond, n = nrow(one), model = NULL, error = conditionMessage(e)))
})

summary_rows <- lapply(fits, function(x) {
  if (is.null(x$model)) return(data.frame(condition = x$condition, n = x$n, Vmax = NA_real_, Vmax_SE = NA_real_, Km = NA_real_, Km_SE = NA_real_, fit_error = x$error))
  estimate <- coef(x$model)
  standard_error <- tryCatch(sqrt(diag(vcov(x$model))), error = function(e) c(Vmax = NA_real_, Km = NA_real_))
  data.frame(condition = x$condition, n = x$n, Vmax = unname(estimate['Vmax']), Vmax_SE = unname(standard_error['Vmax']), Km = unname(estimate['Km']), Km_SE = unname(standard_error['Km']), fit_error = '')
})
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, '/fits.csv', row.names = FALSE, na = '')
writeLines(as.character(dropped_total), '/dropped.txt')

# Base R plot: raw points plus one independently fitted curve for each successful condition.
svg('/plot.svg', width = 8.2, height = 5.2, pointsize = 11, bg = 'white')
if (nrow(d) > 0) {
  palette <- c('#237B69', '#D7874A', '#5878A5', '#A15D73', '#818A48', '#8466A5')
  colors <- rep(palette, length.out = length(conditions))
  plot(d$substrate_mM, d$velocity, type = 'n', log = 'x',
    xlab = 'Substrate concentration (mM)', ylab = expression(Velocity~(mu*mol/min)),
    main = 'Michaelis–Menten kinetics', ylim = range(c(0, d$velocity), finite = TRUE),
    xlim = range(d$substrate_mM, finite = TRUE), las = 1, col.axis = '#47524E',
    col.lab = '#26332F', col.main = '#26332F', bty = 'l')
  grid(col = '#E8ECE8', lty = 1)
  for (i in seq_along(conditions)) {
    cond <- conditions[i]
    points <- d[d$condition == cond, , drop = FALSE]
    points(points$substrate_mM, points$velocity, pch = 19, cex = 1.15, col = colors[i])
    fit <- fits[[i]]
    if (!is.null(fit$model)) {
      estimate <- coef(fit$model)
      curve_x <- exp(seq(log(min(points$substrate_mM)), log(max(points$substrate_mM)), length.out = 250))
      lines(curve_x, estimate['Vmax'] * curve_x / (estimate['Km'] + curve_x), lwd = 2.5, col = colors[i])
    }
  }
  legend('bottomright', legend = conditions, col = colors, pch = 19, lty = 1, lwd = 2, bty = 'n', title = 'Condition')
} else {
  plot.new(); text(0.5, 0.5, 'No valid rows to plot')
}
dev.off()
`;

const $ = (selector) => document.querySelector(selector);
const status = $('#runtime-status');
const runButton = $('#run-button');
const fileInput = $('#csv-file');
const dropZone = $('#drop-zone');
const examplePath = '/data/enzyme_kinetics.csv';
let webR;
let sourceText = '';
let sourceName = 'enzyme_kinetics.csv';
let running = false;

function setStatus(text, ready = false) {
  status.textContent = text;
  status.parentElement.classList.toggle('ready', ready);
}

function showError(message) {
  const box = $('#error-box');
  box.textContent = message;
  box.hidden = false;
}

function hideMessages() {
  $('#error-box').hidden = true;
  $('#warning-box').hidden = true;
}

async function loadExample() {
  setStatus('Loading example data…');
  try {
    const response = await fetch(examplePath);
    if (!response.ok) throw new Error(`Could not load ${examplePath} (${response.status}).`);
    sourceText = await response.text();
    sourceName = 'enzyme_kinetics.csv';
    updateFileInfo();
    $('#results-subtitle').textContent = 'Example dataset ready · separate nls() fits';
    if (webR) runButton.disabled = false;
  } catch (error) {
    setStatus('Example data unavailable');
    showError(`${error.message} You can still choose a CSV file above.`);
  }
}

function updateFileInfo() {
  $('#file-name').textContent = sourceName;
  $('#file-size').textContent = sourceName === 'enzyme_kinetics.csv' ? 'Example dataset' : `${(new Blob([sourceText]).size / 1024).toFixed(1)} KB · ready to analyze`;
}

async function useFile(file) {
  if (!file) return;
  hideMessages();
  try {
    sourceText = await file.text();
    sourceName = file.name;
    updateFileInfo();
    $('#results-subtitle').textContent = 'Custom dataset ready · separate nls() fits';
    runButton.disabled = !webR;
  } catch (error) {
    showError(`Could not read the selected file: ${error.message}`);
  }
}

function decodeFile(file) {
  return new TextDecoder().decode(file);
}

function parseCsv(text) {
  const rows = text.trim().split(/\r?\n/);
  if (rows.length < 2) throw new Error('The CSV must include a header and at least one data row.');
  return rows;
}

function formatEstimate(value, se) {
  if (value === '' || !Number.isFinite(Number(value))) return '<span class="fit-failed">Fit unavailable</span>';
  return `<span class="estimate">${Number(value).toFixed(2)} <small>± ${Number(se).toFixed(2)}</small></span>`;
}

function renderCsvTable(csv, target, type) {
  const rows = parseCsv(csv);
  const headers = rows.shift().split(',');
  const data = rows.filter(Boolean).map((line) => {
    // These output CSVs are generated by R from ordinary condition labels (no quoted commas).
    const values = line.split(',');
    return Object.fromEntries(headers.map((key, index) => [key, values[index] ?? '']));
  });
  if (type === 'counts') {
    target.innerHTML = data.map((row) => `<div class="count-card"><span>${escapeHtml(row.condition)}</span><b>${row.rows_read}<small> rows</small></b><small>${row.rows_used} used</small></div>`).join('');
    return data;
  }
  target.innerHTML = data.map((row) => `<tr><td><span class="condition-label"><i></i>${escapeHtml(row.condition)}</span></td><td>${row.n}</td><td>${formatEstimate(row.Vmax, row.Vmax_SE)}</td><td>${formatEstimate(row.Km, row.Km_SE)}</td></tr>`).join('');
  return data;
}

function escapeHtml(text) {
  return String(text).replace(/[&<>"']/g, (char) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[char]);
}

function buildInterpretation(fits) {
  const box = $('#interpretation');
  const paragraph = box.querySelector('p');
  const control = fits.find((row) => row.condition.toLowerCase() === 'control');
  const inhibitor = fits.find((row) => row.condition.toLowerCase() === 'inhibitor');
  if (!control || !inhibitor || !Number.isFinite(Number(control.Vmax)) || !Number.isFinite(Number(inhibitor.Vmax))) {
    box.hidden = true;
    return;
  }
  const vmaxChange = (Number(inhibitor.Vmax) - Number(control.Vmax)) / Number(control.Vmax);
  const kmChange = (Number(inhibitor.Km) - Number(control.Km)) / Number(control.Km);
  let message = 'The fitted parameter pattern does not clearly match the two simple inhibition patterns listed here.';
  if (Math.abs(vmaxChange) <= 0.2 && kmChange > 0.2) message = 'Vmax is similar while Km increases — consistent with competitive inhibition.';
  else if (vmaxChange < -0.2 && Math.abs(kmChange) <= 0.2) message = 'Vmax decreases while Km is similar — consistent with non-competitive inhibition.';
  paragraph.textContent = message;
  box.hidden = false;
}

async function analyze() {
  if (!webR || !sourceText || running) return;
  running = true;
  runButton.disabled = true;
  runButton.classList.add('loading');
  runButton.querySelector('span:first-child').textContent = 'Running R…';
  hideMessages();
  $('#empty-state').hidden = true;
  $('#results-content').hidden = true;
  $('#results-subtitle').textContent = 'Fitting one nonlinear model per condition…';
  try {
    // Catch missing columns and empty CSVs before sending data to R.
    const header = (sourceText.split(/\r?\n/, 1)[0] || '').split(',').map((item) => item.trim());
    const missing = ['condition', 'substrate_mM', 'velocity'].filter((column) => !header.includes(column));
    if (missing.length) throw new Error(`Missing required column(s): ${missing.join(', ')}.`);
    await webR.FS.writeFile('/input.csv', sourceText);
    await webR.evalRVoid(R_CODE);
    const [fitsRaw, countsRaw, droppedRaw, plotRaw] = await Promise.all([
      webR.FS.readFile('/fits.csv'), webR.FS.readFile('/counts.csv'), webR.FS.readFile('/dropped.txt'), webR.FS.readFile('/plot.svg')
    ]);
    const fitData = renderCsvTable(decodeFile(fitsRaw), $('#fit-table'), 'fits');
    const counts = renderCsvTable(decodeFile(countsRaw), $('#row-counts'), 'counts');
    const dropped = Number(decodeFile(droppedRaw).trim()) || 0;
    const failures = fitData.filter((row) => row.fit_error);
    const warnings = [];
    if (dropped > 0) warnings.push(`${dropped} row${dropped === 1 ? ' was' : 's were'} dropped because condition, substrate, or velocity was missing/invalid (substrate must be > 0).`);
    for (const row of counts) if (Number(row.rows_used) < Number(row.rows_read)) warnings.push(`${row.condition}: ${Number(row.rows_read) - Number(row.rows_used)} row(s) dropped.`);
    for (const row of failures) warnings.push(`${row.condition}: model fit failed — ${row.fit_error}`);
    if (warnings.length) {
      $('#warning-box').innerHTML = warnings.map((item) => `<div>${escapeHtml(item)}</div>`).join('');
      $('#warning-box').hidden = false;
    }
    $('#plot-container').innerHTML = decodeFile(plotRaw);
    $('#r-code').textContent = R_CODE.trim();
    buildInterpretation(fitData);
    $('#results-content').hidden = false;
    $('#results-subtitle').textContent = `${counts.reduce((sum, row) => sum + Number(row.rows_read), 0)} rows read · ${fitData.length} condition${fitData.length === 1 ? '' : 's'} fitted separately`;
  } catch (error) {
    $('#empty-state').hidden = true;
    showError(error.message || String(error));
    $('#results-subtitle').textContent = 'Analysis could not be completed';
  } finally {
    running = false;
    runButton.classList.remove('loading');
    runButton.querySelector('span:first-child').textContent = 'Run analysis';
    runButton.disabled = !sourceText || !webR;
  }
}

runButton.addEventListener('click', analyze);
fileInput.addEventListener('change', (event) => useFile(event.target.files?.[0]));
$('#reset-file').addEventListener('click', () => { hideMessages(); loadExample(); });
dropZone.addEventListener('dragover', (event) => { event.preventDefault(); dropZone.classList.add('dragging'); });
dropZone.addEventListener('dragleave', () => dropZone.classList.remove('dragging'));
dropZone.addEventListener('drop', (event) => { event.preventDefault(); dropZone.classList.remove('dragging'); useFile(event.dataTransfer.files?.[0]); });

async function start() {
  await loadExample();
  try {
    webR = new WebR();
    await webR.init();
    setStatus('Browser R ready', true);
    runButton.disabled = !sourceText;
    $('#results-subtitle').textContent = 'Example dataset ready · separate nls() fits';
  } catch (error) {
    setStatus('Could not start browser R');
    showError(`webR could not initialize: ${error.message || error}. Check your network connection and reload.`);
  }
}

start();
