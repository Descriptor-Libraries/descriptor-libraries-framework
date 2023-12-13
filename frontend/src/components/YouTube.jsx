import React, { useState, useEffect } from 'react';

const ResponsiveYouTube = ({ videoId }) => {
  const [windowWidth, setWindowWidth] = useState(window.innerWidth);

  useEffect(() => {
    const handleResize = () => {
      setWindowWidth(window.innerWidth);
    };

    window.addEventListener('resize', handleResize);

    return () => {
      window.removeEventListener('resize', handleResize);
    };
  }, []);

  // Determine the size of the iframe based on the window width
  const getIframeSize = () => {
    if (windowWidth > 1024) { // Large screens
      return { width: '800px', height: '450px' }; // 16:9 aspect ratio
    } else if (windowWidth > 600) { // Medium screens
      return { width: '600px', height: '338px' }; // 16:9 aspect ratio
    } else { // Small screens
      return { width: '100%', height: 'auto' }; // Full width, auto height
    }
  };

  const iframeSize = getIframeSize();

  return (
    <iframe
      src={`https://www.youtube.com/embed/${videoId}`}
      title="YouTube video player"
      frameBorder="0"
      allowFullScreen
      style={{ width: iframeSize.width, height: iframeSize.height }}
    ></iframe>
  );
};

export default ResponsiveYouTube;
