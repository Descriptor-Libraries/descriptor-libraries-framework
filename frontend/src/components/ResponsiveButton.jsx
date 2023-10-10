import React, { useState, useEffect } from 'react';
import Button from '@mui/material/Button';

function ResponsiveButton({expandedButtonContent, collapsedButtonContent, ...props}) {
  const [width, setWidth] = useState(window.innerWidth);

  useEffect(() => {
    const handleResize = () => {
      setWidth(window.innerWidth);
    };

    window.addEventListener('resize', handleResize);

    return () => {
      window.removeEventListener('resize', handleResize);
    };
  }, []);

  return (
        <Button
          {...props}  
        >
          {width > 768 ? expandedButtonContent : collapsedButtonContent}
        </Button>
      );
}

export default ResponsiveButton;
