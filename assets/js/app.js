/**
 * MIMAS Documentation Website - Main Application JavaScript
 */

// Global state
let flowData = [];
let subroutinesData = {};
let currentStage = 1;

// DOM Elements
let sidebar, stageList, stageSelector, mainContent, searchBox;

// Initialize when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
    initElements();
    loadData();
    setupEventListeners();
});

function initElements() {
    sidebar = document.querySelector('.sidebar');
    stageList = document.querySelector('.stage-list');
    stageSelector = document.querySelector('.stage-selector');
    mainContent = document.querySelector('.main-content');
    searchBox = document.querySelector('.search-box');
}

async function loadData() {
    showLoading();
    
    try {
        // Load flow data
        const flowResponse = await fetch('../data/flow.json');
        if (!flowResponse.ok) throw new Error('Failed to load flow.json');
        flowData = await flowResponse.json();
        
        // Load subroutines data
        const subroutinesResponse = await fetch('../data/subroutines.json');
        if (!subroutinesResponse.ok) throw new Error('Failed to load subroutines.json');
        subroutinesData = await subroutinesResponse.json();
        
        // Render UI
        renderStageList();
        renderStageSelector();
        renderStageContent(currentStage);
        hideLoading();
        
    } catch (error) {
        console.error('Error loading data:', error);
        showCORSError();
    }
}

function setupEventListeners() {
    // Search filter
    if (searchBox) {
        searchBox.addEventListener('input', filterStages);
    }
    
    // Stage selector (mobile)
    if (stageSelector) {
        stageSelector.addEventListener('change', (e) => {
            selectStage(parseInt(e.target.value));
        });
    }
    
    // Accordion click handlers
    document.addEventListener('click', (e) => {
        const header = e.target.closest('.accordion-header');
        if (header) {
            const accordion = header.closest('.accordion');
            toggleAccordion(accordion);
        }
        
        // Snippet button clicks
        if (e.target.classList.contains('snippet-button')) {
            const subName = e.target.dataset.subroutine;
            showSnippet(subName);
        }
    });
    
    // Modal close
    const modal = document.querySelector('.modal-overlay');
    if (modal) {
        modal.addEventListener('click', (e) => {
            if (e.target === modal || e.target.classList.contains('modal-close')) {
                closeModal();
            }
        });
    }
    
    // Keyboard navigation
    document.addEventListener('keydown', (e) => {
        if (e.key === 'Escape') {
            closeModal();
        }
        if (e.key === 'ArrowDown' && !searchBox?.matches(':focus')) {
            e.preventDefault();
            navigateStage(1);
        }
        if (e.key === 'ArrowUp' && !searchBox?.matches(':focus')) {
            e.preventDefault();
            navigateStage(-1);
        }
    });
}

function renderStageList() {
    if (!stageList) return;
    
    stageList.innerHTML = flowData.map(stage => `
        <li class="stage-item" data-stage="${stage.id}" onclick="selectStage(${stage.id})">
            <span class="stage-number">#${stage.id}</span>
            <span class="stage-title">${escapeHtml(stage.title)}</span>
        </li>
    `).join('');
}

function renderStageSelector() {
    if (!stageSelector) return;
    
    stageSelector.innerHTML = flowData.map(stage => `
        <option value="${stage.id}">#${stage.id} - ${escapeHtml(stage.title)}</option>
    `).join('');
}

function renderStageContent(stageId) {
    const stage = flowData.find(s => s.id === stageId);
    if (!stage) return;
    
    // Update active state in sidebar
    document.querySelectorAll('.stage-item').forEach(item => {
        item.classList.toggle('active', parseInt(item.dataset.stage) === stageId);
    });
    
    // Update mobile selector
    if (stageSelector) {
        stageSelector.value = stageId;
    }
    
    // Render content
    mainContent.innerHTML = `
        <div class="content-header">
            <h1>Stage ${stage.id}: ${escapeHtml(stage.title)}</h1>
        </div>
        
        <div class="accordion" id="accordion-a">
            <button class="accordion-header">
                <span class="label">
                    <span class="label-icon a">A</span>
                    What Happens (Overview)
                </span>
                <span class="arrow">▼</span>
            </button>
            <div class="accordion-content">
                <div class="overview-content">
                    ${formatOverview(stage.overview)}
                </div>
            </div>
        </div>
        
        <div class="accordion" id="accordion-b">
            <button class="accordion-header">
                <span class="label">
                    <span class="label-icon b">B</span>
                    Code View
                </span>
                <span class="arrow">▼</span>
            </button>
            <div class="accordion-content">
                <div class="code-content">
                    <div class="section-title">Relevant Subroutines:</div>
                    ${renderSubroutineLinks(stage.relevant_subroutines)}
                    
                    <div class="variable-list">
                        <div class="section-title">Variables Mentioned:</div>
                        ${renderVariables(stage.variables)}
                    </div>
                </div>
            </div>
        </div>
        
        <div class="accordion" id="accordion-c">
            <button class="accordion-header">
                <span class="label">
                    <span class="label-icon c">C</span>
                    Science View
                </span>
                <span class="arrow">▼</span>
            </button>
            <div class="accordion-content">
                <div class="science-content">
                    ${formatScience(stage.science_view)}
                </div>
            </div>
        </div>
    `;
    
    // Open first accordion by default
    document.getElementById('accordion-a')?.classList.add('open');
}

function renderSubroutineLinks(subroutines) {
    if (!subroutines || subroutines.length === 0 || (subroutines.length === 1 && subroutines[0] === 'Not identified yet')) {
        return '<em>No subroutines identified yet</em>';
    }
    
    return subroutines.map(sub => {
        if (sub === 'Not identified yet') return '<em>Not identified yet</em>';
        const hasSnippet = subroutinesData[sub];
        return `
            <span class="subroutine-link" ${hasSnippet ? `data-subroutine="${sub}" onclick="showSnippet('${sub}')"` : ''}>
                ${sub}${hasSnippet ? ' 📄' : ''}
            </span>
        `;
    }).join('');
}

function renderVariables(variables) {
    if (!variables || variables.length === 0 || (variables.length === 1 && variables[0] === 'Not identified yet')) {
        return '<em>No variables identified yet</em>';
    }
    
    return variables.filter(v => v !== 'Not identified yet').map(v => 
        `<span class="variable-item">${v}</span>`
    ).join('');
}

function formatOverview(overview) {
    if (!overview) return '<p>No overview available</p>';
    
    // Convert line breaks to paragraphs
    const paragraphs = overview.split('\n\n').filter(p => p.trim());
    
    return paragraphs.map(p => {
        // Check if it's a list
        if (p.includes('\n')) {
            const lines = p.split('\n').filter(l => l.trim());
            const items = lines.map(l => `<li>${escapeHtml(l.trim())}</li>`).join('');
            return `<ul>${items}</ul>`;
        }
        return `<p>${escapeHtml(p)}</p>`;
    }).join('');
}

function formatScience(science) {
    if (!science) return '<p>No science explanation available</p>';
    
    return `<p>${escapeHtml(science)}</p>`;
}

function filterStages() {
    const query = searchBox.value.toLowerCase();
    
    document.querySelectorAll('.stage-item').forEach(item => {
        const text = item.textContent.toLowerCase();
        item.style.display = text.includes(query) ? '' : 'none';
    });
}

function selectStage(stageId) {
    currentStage = stageId;
    renderStageContent(stageId);
}

function navigateStage(direction) {
    const currentIndex = flowData.findIndex(s => s.id === currentStage);
    const newIndex = currentIndex + direction;
    
    if (newIndex >= 0 && newIndex < flowData.length) {
        selectStage(flowData[newIndex].id);
    }
}

function toggleAccordion(accordion) {
    accordion.classList.toggle('open');
}

function showSnippet(subName) {
    const sub = subroutinesData[subName?.toLowerCase()];
    const modal = document.querySelector('.modal-overlay');
    const modalTitle = document.querySelector('.modal-header h3');
    const modalBody = document.querySelector('.modal-body pre');
    
    if (!modal || !modalTitle || !modalBody) return;
    
    modalTitle.textContent = `SUBROUTINE ${subName}`;
    
    if (sub && sub.snippet) {
        modalBody.textContent = sub.snippet;
    } else {
        modalBody.textContent = `No snippet available for ${subName}`;
    }
    
    modal.classList.add('active');
    document.body.style.overflow = 'hidden';
}

function closeModal() {
    const modal = document.querySelector('.modal-overlay');
    if (modal) {
        modal.classList.remove('active');
        document.body.style.overflow = '';
    }
}

function showLoading() {
    if (mainContent) {
        mainContent.innerHTML = '<div class="loading">Loading data...</div>';
    }
}

function hideLoading() {
    // Loading is replaced by actual content
}

function showCORSError() {
    const message = `
        <div class="error-message">
            <strong>Failed to load data files</strong><br><br>
            This is because browsers block fetch() requests when opening files directly (file:// protocol).<br><br>
            <strong>To fix this, serve the site via HTTP:</strong><br>
            1. Open a terminal in the website folder<br>
            2. Run: <code>python -m http.server 8000</code><br>
            3. Open: <a href="http://localhost:8000" target="_blank">http://localhost:8000</a><br><br>
            <em>The data files (flow.json, subroutines.json) exist and are ready to use.</em>
        </div>
    `;
    if (mainContent) {
        mainContent.innerHTML = message;
    }
    if (stageList) {
        stageList.innerHTML = '<li>' + message + '</li>';
    }
}

function escapeHtml(text) {
    if (!text) return '';
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
}

// Make functions globally available
window.selectStage = selectStage;
window.showSnippet = showSnippet;
