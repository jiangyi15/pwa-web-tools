/**
 * Partial Wave Analysis Tools - Main UI JavaScript
 */

// Run sanity checks when CG module is loaded
document.addEventListener('DOMContentLoaded', () => {
  // Run sanity checks if CG module is available
  if (typeof runSanityChecks === 'function') {
    console.log('');
    runSanityChecks();
  }
  
  // Initialize any interactive components
  initializeAnimations();
});

/**
 * Initialize subtle scroll animations
 */
function initializeAnimations() {
  // Add fade-in animation to cards on scroll
  const observer = new IntersectionObserver((entries) => {
    entries.forEach((entry, index) => {
      if (entry.isIntersecting) {
        setTimeout(() => {
          entry.target.style.opacity = '1';
          entry.target.style.transform = 'translateY(0)';
        }, index * 100);
        observer.unobserve(entry.target);
      }
    });
  }, {
    threshold: 0.1,
    rootMargin: '0px 0px -50px 0px'
  });
  
  document.querySelectorAll('.tool-card').forEach(card => {
    card.style.opacity = '0';
    card.style.transform = 'translateY(20px)';
    card.style.transition = 'opacity 0.5s ease, transform 0.5s ease';
    observer.observe(card);
  });
}

/**
 * Format number with appropriate precision
 */
function formatNumber(num, precision = 10) {
  if (num === 0) return '0';
  
  const abs = Math.abs(num);
  
  if (abs < 0.0001 || abs > 100000) {
    return num.toExponential(precision - 1);
  }
  
  return num.toPrecision(precision);
}

/**
 * Copy text to clipboard
 */
async function copyToClipboard(text) {
  try {
    await navigator.clipboard.writeText(text);
    return true;
  } catch (err) {
    // Fallback for older browsers
    const textArea = document.createElement('textarea');
    textArea.value = text;
    textArea.style.position = 'fixed';
    textArea.style.left = '-9999px';
    document.body.appendChild(textArea);
    textArea.select();
    try {
      document.execCommand('copy');
      return true;
    } catch (e) {
      return false;
    } finally {
      document.body.removeChild(textArea);
    }
  }
}

/**
 * Show a temporary notification
 */
function showNotification(message, type = 'info') {
  const notification = document.createElement('div');
  notification.className = 'notification';
  notification.textContent = message;
  
  const colors = {
    info: 'var(--accent)',
    success: 'var(--success)',
    error: 'var(--error)',
    warning: 'var(--warning)'
  };
  
  notification.style.cssText = `
    position: fixed;
    bottom: 20px;
    right: 20px;
    background: var(--bg-card);
    color: var(--text-primary);
    padding: 0.75rem 1.25rem;
    border-radius: var(--radius-sm);
    border: 1px solid ${colors[type] || colors.info};
    font-size: 0.875rem;
    z-index: 1000;
    animation: slideIn 0.3s ease;
    box-shadow: var(--shadow-lg);
  `;
  
  document.body.appendChild(notification);
  
  setTimeout(() => {
    notification.style.animation = 'slideOut 0.3s ease';
    setTimeout(() => {
      notification.remove();
    }, 300);
  }, 3000);
}

// Add notification animation styles
const style = document.createElement('style');
style.textContent = `
  @keyframes slideIn {
    from {
      opacity: 0;
      transform: translateX(20px);
    }
    to {
      opacity: 1;
      transform: translateX(0);
    }
  }
  
  @keyframes slideOut {
    from {
      opacity: 1;
      transform: translateX(0);
    }
    to {
      opacity: 0;
      transform: translateX(20px);
    }
  }
`;
document.head.appendChild(style);
